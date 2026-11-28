//! Worker threads that compute results for numbered "panels" ahead of the
//! main program (used by RECONR `resxs` and BROADR `broadn`).
//!
//! The main program asks for panel 0, 1, 2, ... in order with `take`. The
//! workers claim panels in the same order and keep at most `window` panels
//! ahead of the one the main program last asked for. A panel's result is
//! only a set of values for the main program to look up; the main program
//! computes anything that is missing itself, so the result of a run never
//! depends on the timing or the number of threads.

use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc::{channel, Receiver};
use std::sync::{Arc, Condvar, Mutex, Once};
use std::thread::JoinHandle;

const THREAD_NAME: &str = "njoy-prefetch";

/// Panics in prefetch threads are not reported (their panel result is just
/// left empty); panics anywhere else are reported as before.
fn quiet_worker_panics() {
    static ONCE: Once = Once::new();
    ONCE.call_once(|| {
        let prev = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            if std::thread::current().name() == Some(THREAD_NAME) {
                return;
            }
            prev(info);
        }));
    });
}

struct Shared {
    /// next panel to claim
    next: AtomicUsize,
    /// panel the main program asked for last
    base: AtomicUsize,
    stop: AtomicBool,
    lock: Mutex<()>,
    cv: Condvar,
}

pub struct RunAhead<R: Send + 'static> {
    sh: Arc<Shared>,
    rx: Receiver<(usize, R)>,
    ready: HashMap<usize, R>,
    handles: Vec<JoinHandle<()>>,
    n: usize,
}

impl<R: Send + Default + 'static> RunAhead<R> {
    /// Start one thread per worker state. `f(worker, p)` computes panel p;
    /// panels first..n are computed.
    pub fn start<W, F>(workers: Vec<W>, first: usize, n: usize, window: usize, f: F) -> RunAhead<R>
    where
        W: Send + 'static,
        F: Fn(&mut W, usize) -> R + Send + Sync + 'static,
    {
        quiet_worker_panics();
        let sh = Arc::new(Shared {
            next: AtomicUsize::new(first),
            base: AtomicUsize::new(first),
            stop: AtomicBool::new(false),
            lock: Mutex::new(()),
            cv: Condvar::new(),
        });
        let f = Arc::new(f);
        let (tx, rx) = channel();
        let mut handles = Vec::new();
        for mut w in workers {
            let sh = sh.clone();
            let f = f.clone();
            let tx = tx.clone();
            let h = std::thread::Builder::new()
                .name(THREAD_NAME.to_string())
                .spawn(move || loop {
                    let q = sh.next.fetch_add(1, Ordering::SeqCst);
                    if q >= n {
                        return;
                    }
                    // stay within the window
                    {
                        let mut g = sh.lock.lock().unwrap();
                        while !sh.stop.load(Ordering::SeqCst) && q > sh.base.load(Ordering::SeqCst) + window {
                            g = sh.cv.wait_timeout(g, std::time::Duration::from_millis(5)).unwrap().0;
                        }
                    }
                    if sh.stop.load(Ordering::SeqCst) {
                        return;
                    }
                    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| f(&mut w, q)));
                    let failed = r.is_err();
                    if tx.send((q, r.unwrap_or_default())).is_err() || failed {
                        // after a panic the worker state may be broken:
                        // this thread stops (its panel was sent empty)
                        return;
                    }
                })
                .expect("prefetch thread");
            handles.push(h);
        }
        RunAhead { sh, rx, ready: HashMap::new(), handles, n }
    }

    /// The result of panel p (waits for it). None if p is out of range or
    /// no worker is left to compute it.
    pub fn take(&mut self, p: usize) -> Option<R> {
        if p >= self.n {
            return None;
        }
        self.sh.base.store(p, Ordering::SeqCst);
        // results of panels before p are no longer needed
        self.ready.retain(|&k, _| k >= p);
        let mut woke = false;
        loop {
            if let Some(r) = self.ready.remove(&p) {
                return Some(r);
            }
            // about to wait: let any worker held back by the window go on
            if !woke {
                let _g = self.sh.lock.lock().unwrap();
                self.sh.cv.notify_all();
                woke = true;
            }
            match self.rx.recv() {
                Ok((q, r)) => {
                    if q >= p {
                        self.ready.insert(q, r);
                    }
                }
                Err(_) => return None,
            }
        }
    }
}

impl<R: Send + 'static> Drop for RunAhead<R> {
    fn drop(&mut self) {
        {
            let _g = self.sh.lock.lock().unwrap();
            self.sh.stop.store(true, Ordering::SeqCst);
            self.sh.cv.notify_all();
        }
        for h in self.handles.drain(..) {
            let _ = h.join();
        }
    }
}
