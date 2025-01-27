// Need to add this dummy interface for the multicore worker to be dependend on
// the external dependencies for mutlicore

use std::marker::PhantomData;

use futures::{
    future::{result, FutureResult},
    Future, IntoFuture, Poll,
};

#[derive(Clone)]
pub struct Worker {
    cpus: usize,
}

impl Worker {
    // We don't expose this outside the library so that
    // all `Worker` instances have the same number of
    // CPUs configured.
    pub(crate) fn new_with_cpus(_cpus: usize) -> Worker {
        Worker { cpus: 1 }
    }

    pub fn new() -> Worker {
        Self::new_with_cpus(1)
    }
    
    #[allow(dead_code)]
    pub fn log_num_cpus(&self) -> u32 {
        0u32
    }

    pub fn compute<F, R>(&self, f: F) -> WorkerFuture<R::Item, R::Error>
    where
        F: FnOnce() -> R + Send + 'static,
        R: IntoFuture + 'static,
        R::Future: Send + 'static,
        R::Item: Send + 'static,
        R::Error: Send + 'static,
    {
        let future = f().into_future();

        WorkerFuture {
            future: result(future.wait()),
        }
    }
    
    #[allow(dead_code)]
    pub fn scope<'a, F, R>(&self, elements: usize, f: F) -> R
    where
        F: FnOnce(&Scope<'a>, usize) -> R,
    {
        let chunk_size = elements;

        let scope = Scope {
            _marker: PhantomData,
        };

        f(&scope, chunk_size)
    }
}
#[derive(Clone)]
pub struct Scope<'a> {
    _marker: PhantomData<&'a usize>,
}

impl<'a> Scope<'a> {
    #[allow(dead_code)]
    pub fn spawn<F, R>(&self, f: F) -> R
    where
        F: FnOnce(&Scope<'a>) -> R,
    {
        f(&self)
    }
}

pub struct WorkerFuture<T, E> {
    future: FutureResult<T, E>,
}

impl<T: Send + 'static, E: Send + 'static> Future for WorkerFuture<T, E> {
    type Item = T;
    type Error = E;

    fn poll(&mut self) -> Poll<Self::Item, Self::Error> {
        self.future.poll()
    }
}
