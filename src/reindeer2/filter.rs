use itertools::Itertools;
use roaring::RoaringBitmap;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
    sync::Mutex,
};

#[cfg(any(debug_assertions, test))]
use std::sync::{atomic::AtomicU64, Arc};

pub struct Filters {
    number_partition: usize,
    data: Vec<Mutex<RoaringBitmap>>,
    color_number: usize,
    bf_bit_size: usize,
    number_abundance: usize,
}

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
            number_partition,
            data,
            color_number,
            bf_bit_size,
            number_abundance,
        }
    }

    pub fn extend_by_draining_partitions_map(
        &self,
        partition_kmers: &mut HashMap<usize, Vec<(u64, u16, usize, usize)>>,
    ) {
        // iterates and empties the hash map when needed
        let partitioned_bf_size = self.bf_bit_size / self.number_partition;
        for (partition_index, kmers) in partition_kmers.drain() {
            let kmer_hashes_to_update =
                kmers
                    .iter()
                    .map(|(kmer_hash, log_abundance, path_num, _chunk_index)| {
                        Self::compute_location(
                            *kmer_hash,
                            partitioned_bf_size,
                            self.color_number,
                            *path_num,
                            self.number_abundance,
                            *log_abundance,
                        )
                    });

            // select the correct BF for the given partition
            let mut bloom_filter = self.data[partition_index]
                .lock()
                .expect("Failed to lock bloom filter");
            // TODO this takes a lot of time
            bloom_filter.extend(kmer_hashes_to_update);
        }
    }

    // TODO should be private
    pub const fn compute_location(
        hash_kmer: u64,
        partitioned_bf_size: usize,
        color_number: usize,      // the total nb of indexed fastas (in a chunk)
        path_color_number: usize, // the index of the current indexed fasta
        abundance_number: usize,
        log_abundance: u16,
    ) -> u32 {
        // compute the position to write
        let location = (hash_kmer % (partitioned_bf_size as u64))
            * (color_number as u64)
            * (abundance_number as u64)
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
    pub fn write_to_disk(
        &self,
        dir_path: &str,
        chunks: &[usize], //  number of colors for each chunk
        partition_nb: usize,
        chunk_id: usize,
    ) -> std::io::Result<()> {
        for (i, bitmap) in self.data.iter().enumerate() {
            let p = i % partition_nb;
            assert_eq!(p, i);
            let file_path = Path::new(dir_path)
                .join(format!("partition_bloom_filters_c{}_p{}.bin", chunk_id, p));
            let file = File::create(&file_path)?;
            let mut writer = BufWriter::new(file);
            let chunk_colors = chunks[chunk_id];
            // write the number of colors as a u64 to the file first
            writer.write_all(&chunk_colors.to_le_bytes())?;
            // serialize the bitmap into the file
            let mut locked_bitmap = bitmap.lock().expect("fatal error: a thread holding the mutex panicked, so this thread will panic as well");
            locked_bitmap.serialize_into(&mut writer)?;
            locked_bitmap.clear();
        }

        Ok(())
    }
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
