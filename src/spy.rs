
pub struct SpyIter<I> where I: Iterator, I::Item: Copy {
    iter: I,
    count: usize,
    done: bool,
    buf: Vec<I::Item>,
    callback: Option<Box<dyn FnOnce(Vec<I::Item>)>>,
}

impl <I> Iterator for SpyIter<I> where I: Iterator, I::Item: Copy {
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        match self.iter.next() {
            None => {
                if !self.done {
                    self.done = true;
                    self.callback.take().map(|cb| cb(std::mem::take(&mut self.buf)));
                }
                None
            }
            Some(item) => {
                if !self.done {
                    self.buf.push(item);
                    if self.buf.len() >= self.count {
                        self.done = true;
                        self.callback.take().map(|cb| cb(std::mem::take(&mut self.buf)));
                    }
                }
                Some(item)
            },
        }
    }
}

pub trait SpyableIter: Iterator + Sized {
    fn spy<T>(self, count: usize, callback: T) -> SpyIter<Self> where T: FnOnce(Vec<Self::Item>) + 'static, Self::Item: Copy {
        SpyIter {
            iter: self,
            count,
            done: false,
            buf: Vec::with_capacity(count),
            callback: Some(Box::new(callback)),
        }
    }
}

impl<I> SpyableIter for I where I: Iterator + Sized {}
