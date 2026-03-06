// deny the use of unwrap: use expect instead
#![cfg_attr(not(test), deny(clippy::unwrap_used))]
#![deny(dead_code)]

pub mod reindeer2;
