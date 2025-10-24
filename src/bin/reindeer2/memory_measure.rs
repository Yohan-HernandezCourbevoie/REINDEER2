use libc::{getrusage, rusage, RUSAGE_SELF};
use std::mem::MaybeUninit;

pub fn format_int_with_spaces(mut n: i64) -> String {
    if n == 0 {
        return "0".to_string();
    }
    let negative = n < 0;
    if negative {
        n = -n;
    }

    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    let mut cnt = 0;
    for ch in s.chars().rev() {
        if cnt == 3 {
            out.push(' ');
            cnt = 0;
        }
        out.push(ch);
        cnt += 1;
    }
    if negative {
        out.push('-');
    }
    out.chars().rev().collect()
}

pub fn get_max_rss() -> i64 {
    let mut usage = MaybeUninit::<rusage>::uninit();
    let usage = unsafe {
        getrusage(RUSAGE_SELF, usage.as_mut_ptr());
        usage.assume_init()
    };
    usage.ru_maxrss
}
