use std::{
    fs::File,
    io::{self, BufWriter, Read, Write},
    path::Path,
    sync::{Arc, Mutex},
};

pub fn write_kmer_counts_to_disk(
    dir_path: &Path,
    kmer_counts_vector: &Arc<Mutex<Vec<usize>>>,
) -> io::Result<()> {
    // OPTIMIZE we may be able to drop the lock before writing to disk
    let file_path = dir_path.join("kmer_counts_per_color.bin");
    let file = File::create(&file_path)?;
    let mut writer = BufWriter::new(file);
    let mut locked_vector = kmer_counts_vector.lock().expect(
        "fatal error: a thread holding the mutex panicked, so this thread will panic as well",
    );
    let locked_vector_ref: &Vec<usize> = &locked_vector;
    let binary_encoded = bincode::serialize(locked_vector_ref)
        .expect("should have been able to serialize the count of k-mers");
    locked_vector.clear();
    drop(locked_vector);
    writer.write_all(&binary_encoded)?;
    Ok(())
}

pub fn load_kmer_counts_vector(dir_path: &Path) -> io::Result<Vec<usize>> {
    let mut file = File::open(dir_path.join("kmer_counts_per_color.bin"))?;
    // Read the rest of the file to deserialize the hashmap
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    let counts_vector = bincode::deserialize_from(&buffer[..]).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            "Failed to deserialize the counts vector",
        )
    })?;
    Ok(counts_vector)
}
