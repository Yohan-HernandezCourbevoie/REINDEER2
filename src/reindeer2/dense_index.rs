use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Read, Write},
    path::Path,
    sync::{
        atomic::{AtomicU64, Ordering},
        Mutex,
    },
};

use crate::reindeer2::filter::Filters;
use serde::{Deserialize, Serialize};

fn count_zeros(abundance_vector: &[u8], max_index: usize) -> usize {
    abundance_vector
        .iter()
        .take(max_index + 1)
        .filter(|&&val| val == 0)
        .count()
}

#[derive(Serialize, Deserialize)]
pub struct DenseIndexPartition {
    data: HashMap<u64, Vec<u8>>,
}

impl DenseIndexPartition {
    pub fn new() -> Self {
        Self {
            data: HashMap::new(),
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            data: HashMap::with_capacity(capacity),
        }
    }

    pub fn insert(&mut self, kmer_hash: u64, abundance_vector: Vec<u8>) -> Option<Vec<u8>> {
        self.data.insert(kmer_hash, abundance_vector)
    }

    pub fn contains_key(&self, kmer_hash: &u64) -> bool {
        self.data.contains_key(kmer_hash)
    }

    pub fn get_mut(&mut self, kmer_hash: &u64) -> Option<&mut Vec<u8>> {
        self.data.get_mut(kmer_hash)
    }

    pub fn dump(&mut self, writer: &mut impl Write) -> std::io::Result<()> {
        let binary_encoded = bincode::serialize(self).unwrap();
        writer.write_all(&binary_encoded)?;
        self.data.clear();
        Ok(())
    }

    pub fn get_abundance(&self, kmer: &u64) -> Option<&Vec<u8>> {
        self.data.get(kmer)
    }

    pub fn load_from_disk(file_path: &str) -> std::io::Result<Self> {
        let mut file = File::open(file_path)?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        bincode::deserialize_from(&buffer[..]).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Failed to deserialize hashmap",
            )
        })
    }

    // TODO name
    pub fn muche(
        &mut self,
        bloom_filters: &Filters,
        path_num: usize,
        threshold: usize,
        max_map_size: usize,
        atomic_dense_kmers_count: &AtomicU64,
        atomic_sparse_kmers_count: &AtomicU64,
        chunk_index: usize,
        partition_kmers: &mut HashMap<usize, Vec<(u64, u16, usize, usize)>>,
        partition_index: usize,
    ) {
        let mut kmers_to_remove: Vec<u64> = Vec::new();
        for (kmer_hash, abundance_vector) in self.data.iter() {
            let number_of_zeros = count_zeros(abundance_vector, path_num);
            if number_of_zeros > threshold {
                let color_count = path_num - number_of_zeros + 1;
                atomic_dense_kmers_count.fetch_sub(color_count as u64, Ordering::Relaxed);
                atomic_sparse_kmers_count.fetch_add(color_count as u64, Ordering::Relaxed);
                kmers_to_remove.push(*kmer_hash);
                for (path_index, log_abundance) in
                    abundance_vector.iter().take(path_num).enumerate()
                {
                    if *log_abundance > 0 {
                        partition_kmers // separate the kmers per partition
                            .entry(partition_index)
                            .or_default()
                            .push((
                                *kmer_hash,
                                (*log_abundance - 1) as u16,
                                path_index,
                                chunk_index,
                            ));
                    }
                }
                if partition_kmers.len() >= max_map_size {
                    bloom_filters.extend_by_draining_partitions_map(partition_kmers);
                }
            }
        }
        kmers_to_remove.iter().for_each(|kmer| {
            self.data.remove(kmer);
        });
        self.data.shrink_to_fit();
    }
}

#[derive(Serialize, Deserialize)]
pub struct DenseIndex {
    data: Vec<Mutex<DenseIndexPartition>>,
}

impl DenseIndex {
    pub fn with_partition_number(partition_number: usize) -> Self {
        let data = (0..partition_number)
            // TODO magic number here
            .map(|_| {
                Mutex::new(DenseIndexPartition::with_capacity(
                    200_000_000 / partition_number,
                ))
            })
            .collect();
        Self { data }
    }

    pub fn write_to_disk(&self, dir_path: &str) -> std::io::Result<()> {
        for (partition_id, partition) in self.data.iter().enumerate() {
            let file_path =
                Path::new(dir_path).join(format!("partition_dense_index_p{}.bin", partition_id));
            let file = File::create(&file_path)?;
            let mut writer: BufWriter<File> = BufWriter::new(file);
            let mut partition = partition.lock().unwrap();
            partition.dump(&mut writer)?;
        }
        Ok(())
    }

    // pub fn insert_abundance_to_partition(
    //     &self,
    //     partition_index: usize,
    //     kmer_hash: u64,
    //     abundances: Vec<u8>,
    // ) {
    //     self.data[partition_index]
    //         .lock()
    //         .expect("Failed to lock the dense index")
    //         .insert(kmer_hash, abundances);
    // }

    /// Inserts in the index if the k-mer can be dense.
    /// Returns `true` if the k-mer was inserted, `false` otherwise
    pub fn insert_if_dense(
        &self,
        partition_index: usize,
        kmer_hash: u64,
        path_num_global: usize,
        threshold: usize,
        log_abundance: u16,
        color_number_global: usize,
    ) -> bool {
        let mut partition = self.data[partition_index]
            .lock()
            .expect("Failed to lock the dense index");
        if partition.contains_key(&kmer_hash) {
            // update the vector with the right abundance
            if let Some(abundance_vector) = partition.get_mut(&kmer_hash) {
                abundance_vector[path_num_global] = (log_abundance + 1) as u8;
            }
            return true;
        }
        if path_num_global <= threshold {
            // create a new abundance vector for the k-mer
            let mut abundance_vector: Vec<u8> = vec![0; color_number_global];
            abundance_vector[path_num_global] = (log_abundance + 1) as u8;
            partition.insert(kmer_hash, abundance_vector);
            return true;
        }
        false
    }

    // TODO name
    pub fn truc(
        &self,
        bloom_filters: &Filters,
        path_num: usize,
        threshold: usize,
        max_map_size: usize,
        atomic_dense_kmers_count: &AtomicU64,
        atomic_sparse_kmers_count: &AtomicU64,
        chunk_index: usize,
    ) {
        let mut partition_kmers: HashMap<usize, Vec<(u64, u16, usize, usize)>> = HashMap::new(); // keep kmer info to fill BFs
        for (partition_index, hashmap) in self.data.iter().enumerate() {
            let mut dense_index = hashmap // select the correct BF for the given partition
                .lock()
                .expect("Failed to lock bloom filter");
            dense_index.muche(
                bloom_filters,
                path_num,
                threshold,
                max_map_size,
                atomic_dense_kmers_count,
                atomic_sparse_kmers_count,
                chunk_index,
                &mut partition_kmers,
                partition_index,
            );
        }
        // Flush remaining k-mers in the map by calling the earlier closure
        bloom_filters.extend_by_draining_partitions_map(&mut partition_kmers);
    }
}
