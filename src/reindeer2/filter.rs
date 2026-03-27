use itertools::Itertools;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use std::{
    cmp::min,
    collections::HashMap,
    fs::File,
    io::{BufReader, BufWriter, Read, Write},
    path::Path,
    sync::Mutex,
};

use crate::reindeer2::{create_and_reserve_tar_get_file, NB_FILE_IN_AN_INDEX};

#[cfg(any(debug_assertions, test))]
use std::sync::{atomic::AtomicU64, Arc};

pub struct Filters {
    data: Vec<Mutex<RoaringBitmap>>,
    color_number: usize,
    bf_bit_size: usize,
    number_abundance: usize,
}

const LOW_BITS: u64 = 16;
const MINIMIZER_BIT_MASK: u64 = u64::MAX << LOW_BITS;
const KMER_BIT_MASK: u64 = !MINIMIZER_BIT_MASK;
// const BLOCK_SIZE_BIT: u64 = 1 << (LOW_BITS - 1);

impl Filters {
    pub fn with_number_partition(
        number_partition: usize,
        color_number: usize,
        bf_bit_size: usize,
        number_abundance: usize,
    ) -> Self {
        let new_mutex = |_| Mutex::new(RoaringBitmap::new());
        let data = (0..number_partition).map(new_mutex).collect_vec();

        Self {
            data,
            color_number,
            bf_bit_size,
            number_abundance,
        }
    }

    pub fn extend_by_draining_partitions_map(
        &self,
        partition_kmers: &mut HashMap<usize, Vec<(u64, u64, u16, usize, usize)>>,
    ) {
        // iterates and empties the hash map when needed
        let partitioned_bf_size = self.bf_bit_size / self.data.len();
        for (partition_index, kmers) in partition_kmers.drain() {
            let kmer_hashes_to_update = kmers.iter().map(
                |(kmer_hash, minimizer_hash, log_abundance, path_num, _chunk_index)| {
                    Self::compute_location(
                        *kmer_hash,
                        *minimizer_hash,
                        partitioned_bf_size,
                        self.color_number,
                        *path_num,
                        self.number_abundance,
                        *log_abundance,
                    )
                },
            );

            // select the correct BF for the given partition
            let mut bloom_filter = self.data[partition_index]
                .lock()
                .expect("Failed to lock bloom filter");
            // TODO this takes a lot of time
            bloom_filter.extend(kmer_hashes_to_update);
        }
    }

    /// Returns the first position of the concerned column in the partition
    pub const fn compute_base_position(
        smer_hash: u64,
        minimizer_hash: u64,
        partitioned_bf_size: usize,
        color_number: usize,
        abundance_number: usize,
    ) -> u64 {
        let nb_block = partitioned_bf_size as u64 / LOW_BITS;
        let block_id = minimizer_hash % nb_block;
        let block_start = block_id * LOW_BITS;
        let smer_hash = (block_start & MINIMIZER_BIT_MASK) | (smer_hash & KMER_BIT_MASK);

        let position = smer_hash % (partitioned_bf_size as u64);
        position * (color_number as u64) * (abundance_number as u64)
    }

    // TODO should be private
    pub const fn compute_location(
        smer_hash: u64,
        minimizer_hash: u64,
        partitioned_bf_size: usize,
        color_number: usize,      // the total nb of indexed fastas (in a chunk)
        path_color_number: usize, // the index of the current indexed fasta
        abundance_number: usize,
        log_abundance: u16,
    ) -> u32 {
        // compute the position to write
        let base_position = Self::compute_base_position(
            smer_hash,
            minimizer_hash,
            partitioned_bf_size,
            color_number,
            abundance_number,
        );
        let location = base_position
            + (path_color_number as u64) * (abundance_number as u64)
            + (log_abundance as u64);
        debug_assert!(location == (location as u32) as u64);
        location as u32

        /* example
                            c0  c1  c2  c3
            color 0   abund 0   0   1   1
                      abund 1   0   0   0
            color 1   abund 0   0   1   0
                      abund 1   0   0   1
            color 2   abund 0   0   0   1
                      abund 1   0   1   0
            with two partitions becomes two files
            f0=010101000000

            f1=101001100110

            let's say a k-mer x is in color 2 with abund 0, hashed in column 2
                            c0  c1  c2  c3
            color 0   abund 0   0   1   1
                      abund 1   0   0   0
            color 1   abund 0   0   1   0
                      abund 1   0   0   1
            color 2   abund 0   0   X   1
                      abund 1   0   1   0

            so we must do a modification so f1 becomes
            f1=1=1010X1100110 with X = 1, so index 4 in the vector

            partition number is 1, and the partition size (number of columns) is 2, composed of column 2 and 3.
            We'll have hash(x)%2=0, which means we're in the first chunk of the vector ((hash_k % partitioned_bf_size)*color_number*abundance_number = 0 here)
            then we want to go up to color 2, so we mus pass through color 0 and 1 that have two abundances each (path_color_number*abundance_number = 2*2 here)
            then go to the right abundance, here 0 that corresponds to log_abundance
            position_to_write = 0 + 2*2 + 0 = index 4 in the vector
        */
    }

    // TODO take a Write (or even better, implement a trait (Serialize ?)
    pub fn write_to_disk_as_multiple_small_files(
        &self,
        dir_path: &str,
        chunks: &[usize], //  number of colors for each chunk
        chunk_id: usize,
    ) -> std::io::Result<()> {
        let chunk_colors = chunks[chunk_id];
        for (i, bitmap) in self.data.iter().enumerate() {
            let file_path = Path::new(dir_path)
                .join(format!("partition_bloom_filters_c{}_p{}.bin", chunk_id, i));
            let file = File::create(&file_path)?;
            let mut writer = BufWriter::new(file);
            // write the number of colors as a u64 to the file first
            writer.write_all(&chunk_colors.to_le_bytes())?;
            // serialize the bitmap into the file
            let mut locked_bitmap = bitmap.lock().expect("fatal error: a thread holding the mutex panicked, so this thread will panic as well");
            locked_bitmap.serialize_into(&mut writer)?;
            locked_bitmap.clear();
        }

        Ok(())
    }

    pub fn write_to_disk_tar_get_files(self, output_dir: &str) -> std::io::Result<()> {
        let nb_partition = self.data.len();
        let nb_partition_in_a_file = nb_partition.div_ceil(NB_FILE_IN_AN_INDEX);
        (0..NB_FILE_IN_AN_INDEX)
            .into_par_iter() // TODO
            .for_each(|file_id| {
                let range_start = nb_partition_in_a_file * file_id;
                let range_end = min(nb_partition_in_a_file * (file_id + 1), nb_partition);
                let range = range_start..range_end;
                let path = Path::new(output_dir)
                    .join(format!("partition_bloom_filters_group{file_id}.bin"));
                let len: usize = range
                    .try_len()
                    .expect("range object should have a known length");
                let len = u64::try_from(len)
                    .expect("should have less than u64::MAX partitions in a single file");
                let mut file = create_and_reserve_tar_get_file(&path, len);
                for i in range {
                    let bitmap = &self.data[i];
                    let mut bitmap = bitmap.lock().expect("fatal error: a thread holding the mutex panicked, so this thread will panic as well");
                    tar_get::append_element(&mut file, &bitmap, |writer, bitmap| {
                        bitmap.serialize_into(writer)
                    })
                    .expect("should have been able to add an element to a tar_get file");
                    bitmap.clear();
                }
            });
        Ok(())
    }
}

/// Loads a Bloom filter from disk.
pub fn load_bloom_filter_from_big_file(
    file_path: &str,
    index: u64,
) -> std::io::Result<RoaringBitmap> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let deserializer = |reader: Vec<u8>| {
        let slice: &[u8] = &reader;
        RoaringBitmap::deserialize_from(slice)
    };
    let bitmap = tar_get::deserialize(reader, index, deserializer)
        .unwrap_or_else(|_| panic!("should have been able to load index {index} from {file_path}")); // TODO error handling
    Ok(bitmap)
}

/// Loads a raw Bloom filter from disk
pub fn load_raw_bloom_filter(file_path: &str) -> std::io::Result<(RoaringBitmap, usize)> {
    let file = File::open(file_path)?;
    let mut reader = BufReader::new(file);

    // read the first 8 bytes as a u64 to get the number of colors
    let mut color_buffer = [0u8; 8];
    reader.read_exact(&mut color_buffer)?;
    let local_color_nb = u64::from_le_bytes(color_buffer) as usize;

    // Rread the rest of the file to deserialize the Bloom filter
    let mut buffer = Vec::new();
    reader.read_to_end(&mut buffer)?;
    let bitmap = RoaringBitmap::deserialize_from(&buffer[..]).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Failed to deserialize bitmap",
        )
    })?;
    Ok((bitmap, local_color_nb))
}

#[cfg(any(debug_assertions, test))]
impl Filters {
    // TODO document this function
    // TODO this are outparameters: explain in the docs
    pub fn update_sparse_counts(
        &self,
        atomic_sparse_one_seen: &AtomicU64,
        atomic_sparse_fp_seen: &AtomicU64,
        abundance_number: usize,
    ) -> std::io::Result<()> {
        use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
        use std::sync::atomic::Ordering;
        // use std::{
        //     fs::File,
        //     io::{BufWriter, Write},
        // };

        let ones_by_partition = Arc::new(Mutex::new(Vec::<(usize, usize)>::new()));
        self.data.par_iter().for_each(|bitmap| {
            let locked_bitmap = bitmap.lock().expect("fatal error: a thread holding the mutex panicked, so this thread will panic as well");
            let mut color: usize = 0;
            let mut new_color: usize = 0;
            let mut ones: usize = 0;
            let mut fp: usize = 0;
            locked_bitmap.iter().for_each(|value| {
                ones += 1;
                new_color = (value as usize) / abundance_number;
                if new_color != color {
                    color = new_color;
                } else {
                    fp += 1;
                }
            });
            let mut tmp_vector = ones_by_partition.lock().expect("fatal error: a thread holding the mutex panicked, so this thread will panic as well");
            tmp_vector.push((ones, fp));
            atomic_sparse_one_seen.fetch_add(ones as u64, Ordering::Relaxed);
            atomic_sparse_fp_seen.fetch_add(fp as u64, Ordering::Relaxed);
        });

        // let file = File::create("partitions_info.log")?;
        // let mut writer = BufWriter::new(file);
        // let tmp_vector: std::sync::MutexGuard<'_, Vec<(usize, usize)>> = ones_by_partition
        //     .lock()
        //     .expect(
        //     "fatal error: a thread holding the mutex panicked, so this thread will panic as well",
        // );
        // for (a1, a2) in tmp_vector.iter() {
        //     writer.write_all(format!("{},{}\n", a1, a2).as_bytes())?;
        // }
        // writer.flush()?;

        Ok(())
    }
}
