pub struct Interleaver<I,T> where I: Iterator<Item = (T,T)> {
    iter: I,
    is_l: bool,
    r: T
}

impl<I,T> Interleaver<I,T> where I: Iterator<Item = (T,T)>, T: Copy + Default {
    fn new(iter: I) -> Self {
        Interleaver {
            iter,
            is_l: true,
            r: T::default()
        }
    }
}

impl<I,T> Iterator for Interleaver<I,T> where I: Iterator<Item = (T,T)>, T: Copy + Default {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_l {
            let (l, r) = self.iter.next()?;
            self.r = r;
            self.is_l = false;
            Some(l)
        } else {
            self.is_l = true;
            Some(self.r)
        }
    }
}

pub trait InterleaveableIter<T> where T: Copy + Default {
    fn interleave(self) -> Interleaver<Self, T> where Self: Sized, Self: Iterator<Item = (T,T)> {
        Interleaver::new(self)
    }
}

impl<I,T> InterleaveableIter<T> for I where I: Iterator<Item = (T,T)> + Sized, T: Copy + Default {}