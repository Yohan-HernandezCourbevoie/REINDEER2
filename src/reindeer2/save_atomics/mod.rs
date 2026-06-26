pub mod atomics;

#[cfg(any(debug_assertions, test))]
pub mod debug_atomics;
