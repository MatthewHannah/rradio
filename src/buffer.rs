use std::{
    ops::{Deref, DerefMut},
    sync::mpsc,
    fmt::Debug
};

enum QueueItem<T> {
    Item(T),
    Done
}

#[derive(Debug)]
pub struct BufToken<T>
where
    T: Sized + Default,
{
    data: T,
    released: bool,
    ret: mpsc::Sender<QueueItem<T>>,
}

impl<T> Deref for BufToken<T>
where
    T: Default,
{
    type Target = T;
    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T> DerefMut for BufToken<T>
where
    T: Default,
{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

impl<T> Drop for BufToken<T>
where
    T: Default,
{
    fn drop(&mut self) {
        if !self.released {
            // ignore errors if we need to here
            let _ = self.ret.send(QueueItem::Item(std::mem::take(&mut self.data)));
        }
        self.released = true;
    }
}

#[derive(Debug)]
pub struct SendBuf<T>
where
    T: Default,
{
    buf_in: mpsc::Receiver<QueueItem<T>>,
    buf_out: mpsc::Sender<QueueItem<T>>,
    buf_return: mpsc::Sender<QueueItem<T>>,
    done: bool, // used since we have a circular reference to 
}

impl<T> SendBuf<T>
where
    T: Default,
{
    pub fn get(&mut self) -> Option<BufToken<T>> {
        // if we know we are done, just return
        if self.done {
            None
        } else {
            // try to get a T from the pool
            // should always return Ok() bercause we have a circular reference with buf_in/buf_return
            // but we check anyway
            if let Ok(t) = self.buf_in.recv() {
                match t {
                    QueueItem::Item(t) => Some(BufToken {
                        data: t,
                        released: false,
                        ret: self.buf_return.clone(),
                    }),
                    QueueItem::Done => {
                        self.done = true;
                        None
                    },
                }
            } else {
                None
            }
        }
    }

    pub fn commit(&mut self, mut token: BufToken<T>) {
        let _ = self.buf_out.send(QueueItem::Item(std::mem::take(&mut token.data)));
        token.released = true;
    }
}

impl<T> Drop for SendBuf<T> where T: Default {
    fn drop(&mut self) {
        let _ = self.buf_out.send(QueueItem::Done);
    }
}

#[derive(Debug)]
pub struct RecvBuf<T>
where
    T: Default,
{
    buf_in: mpsc::Receiver<QueueItem<T>>,
    buf_out: mpsc::Sender<QueueItem<T>>,
}

impl<T> RecvBuf<T>
where
    T: Default,
{
    pub fn get(&mut self) -> Option<BufToken<T>> {
        match self.buf_in.recv() {
            Ok(QueueItem::Item(t)) => Some(BufToken {
                data: t,
                released: false,
                ret: self.buf_out.clone(),
            }),
            Ok(QueueItem::Done) => None,
            Err(_) => None
        }
    }

    pub fn release(&mut self, mut token: BufToken<T>) {
        let _ = self.buf_out.send(QueueItem::Item(std::mem::take(&mut token.data)));
        token.released = true;
    }
}

impl<T> Drop for RecvBuf<T> where T: Default {
    fn drop(&mut self) {
        let _ = self.buf_out.send(QueueItem::Done);
    }
}

pub fn buf_pair<T>(capacity: usize) -> (SendBuf<T>, RecvBuf<T>)
where
    T: Default,
{
    let (buf_ret_tx, buf_ret_rx) = mpsc::channel();
    let (buf_tx, buf_rx) = mpsc::channel();

    for _ in 0..capacity {
        buf_ret_tx.send(QueueItem::Item(T::default())).unwrap();
    }

    let send = SendBuf {
        buf_in: buf_ret_rx,
        buf_out: buf_tx,
        buf_return: buf_ret_tx.clone(),
        done: false,
    };

    let recv = RecvBuf {
        buf_in: buf_rx,
        buf_out: buf_ret_tx,
    };

    (send, recv)
}

#[derive(Debug)]
pub struct RecvBufIter<T> where T: Copy {
    rx: RecvBuf<Vec<T>>,
    curr: Option<BufToken<Vec<T>>>,
    idx: usize,
}

impl<T> RecvBufIter<T> where T: Copy {
    pub fn new(buf: RecvBuf<Vec<T>>) -> RecvBufIter<T> {
        RecvBufIter { rx: buf, curr: None, idx: 0 }
    }

    fn maybe_load_next(&mut self) -> Option<()> {
        // below keeps us looping if we get passed zero length buffer
        while self.curr.is_none() || self.idx >= self.curr.as_ref()?.len() {
            self.idx = 0;

            if self.curr.is_some() {
                self.rx.release(self.curr.take().unwrap()); // we are done with this buffer
            }

            // if RecvBuf says we are done, return None
            if let Some(tok) = self.rx.get() {
                self.curr = Some(tok);
            } else {
                return None;
            }
        }

        Some(())
    }
}

impl<T> Iterator for RecvBufIter<T> where T: Copy + Debug {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr.is_none() || self.idx >= self.curr.as_ref().expect("should be some, checked is none").len() {
            self.maybe_load_next()?;
        }
        
        let idx: usize = self.idx;    
        self.idx += 1;
        Some(self.curr.as_ref()?[idx])
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_buf() {
        let (mut tx, mut rx) = buf_pair::<u32>(5);

        let mut token = tx.get().unwrap();
        *token = 3;
        tx.commit(token);

        let token = rx.get().unwrap();
        assert_eq!(*token, 3);
        rx.release(token);

        for i in 0..5 {
            let mut token = tx.get().unwrap();
            *token = i;
            tx.commit(token);
        }

        for i in 0..5 {
            let token = rx.get().unwrap();
            assert_eq!(i, *token);
            rx.release(token);
        }
    }

    #[test]
    fn test_rx_drop() {
        let (mut tx, mut rx) = buf_pair::<u32>(5);

        let mut token = tx.get().unwrap();
        *token = 3;
        tx.commit(token);

        let token = rx.get().unwrap();
        assert_eq!(*token, 3);
        rx.release(token);

        for i in 0..5 {
            let mut token = tx.get().unwrap();
            *token = i;
            tx.commit(token);
        }

        for i in 0..5 {
            let token = rx.get().unwrap();
            assert_eq!(i, *token);

            // intentionally do not release, should still return the buf after token drops
            // rx.release(token);
        }
        
        // should be able to still reacquire buffers to send now
        for i in 0..5 {
            let mut token = tx.get().unwrap();
            *token = i;
            tx.commit(token);
        }
    }

    #[test]
    fn test_tx_drop() {
        let (mut tx, mut rx) = buf_pair::<u32>(5);

        let mut token = tx.get().unwrap();
        *token = 3;
        tx.commit(token);

        let token = rx.get().unwrap();
        assert_eq!(*token, 3);
        rx.release(token);

        for i in 0..4 {
            let mut token = tx.get().unwrap();
            *token = i;
            tx.commit(token);
        }

        {
            let mut token = tx.get().unwrap();
            *token = 5;
        }

        {
            let mut token = tx.get().unwrap();
            *token = 6;
            tx.commit(token);
        }

        for i in 0..4 {
            let token = rx.get().unwrap();
            assert_eq!(i, *token);
            rx.release(token);
        }

        let token = rx.get().unwrap();
        assert_eq!(6, *token);
    }

    use std::thread;
    #[test]
    fn test_basic_threading() {
        let (mut tx, mut rx) = buf_pair::<Vec<u32>>(2);

        let rx_handle = thread::spawn(move || {
            for i in 0..50 {
                let token = rx.get().unwrap();
                assert!(token.len() == i);
                rx.release(token);
            }
        });

        let tx_handle = thread::spawn(move || {
            for i in 0..50 {
                let mut token = tx.get().unwrap();
                token.resize(i, i as u32);
                tx.commit(token);
            }
        });

        tx_handle.join().unwrap();
        rx_handle.join().unwrap();
    }

    #[test]
    fn test_threading_rx_early_end() {
        let (mut tx, mut rx) = buf_pair::<Vec<u32>>(2);

        let rx_handle = thread::spawn(move || {
            for i in 0..10 {
                let token = rx.get().unwrap();
                assert!(token.len() == i);
                rx.release(token);
            }
        });

        let tx_handle = thread::spawn(move || {
            for i in 0..50 {
                if let Some(mut token) = tx.get() {
                    token.resize(i, i as u32);
                    tx.commit(token);
                } else {
                    break;
                }
            }
        });

        tx_handle.join().unwrap();
        rx_handle.join().unwrap();
    }

    #[test]
    fn test_threading_tx_early_end() {
        let (mut tx, mut rx) = buf_pair::<Vec<u32>>(2);

        let rx_handle = thread::spawn(move || {
            for i in 0..50 {
                if let Some(token) = rx.get() {
                    assert!(token.len() == i);
                    rx.release(token);
                } else {
                    break;
                }
            }
        });

        let tx_handle = thread::spawn(move || {
            for i in 0..10 {
                let mut token = tx.get().unwrap();
                token.resize(i, i as u32);
                tx.commit(token);
            }
        });

        tx_handle.join().unwrap();
        rx_handle.join().unwrap();
    }

    #[test]
    fn test_rx_iter() {
        let (mut tx, rx) = buf_pair::<Vec<u32>>(2);

        let rx_handle = thread::spawn(move || {
            let rx_iter = RecvBufIter::new(rx);
            let mut any = 0;
            for (c, x) in rx_iter.enumerate() {
                assert!(x == (c as u32));
                any += 1;
            }

            assert_ne!(any, 0);
        });

        let tx_handle = thread::spawn(move || {
            let mut counter = 0;
            for i in 0..10 {
                if let Some(mut token) = tx.get() {
                    token.resize(i, i as u32);
                    for c in token.iter_mut() {
                        *c = counter;
                        counter += 1;
                    }
                    tx.commit(token);
                } else {
                    break;
                }
            }
        });

        tx_handle.join().unwrap();
        rx_handle.join().unwrap();
    }

}
