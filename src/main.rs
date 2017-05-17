extern crate rand;
extern crate advapi32;
extern crate kernel32;
extern crate winapi;
extern crate time;
use rand::{StdRng, SeedableRng, Rng};
use std::ops::{Index, IndexMut};
use std::cmp::Ordering;
use std::io::Write;
use time::precise_time_ns;

const RAND_MAX: u32 = 32768;

const N: usize = 1200;


fn main() {
    let tstart = precise_time_ns();
    let arena = large_page_alloc::init();

    let mut rng: StdRng = SeedableRng::from_seed(&[0usize] as &[usize]);

    // arena.alloc returns reference to uninitialized value
    let mut a = unsafe{ arena.alloc::<[f64;N*N]>() };
    let mut a = A2DSlice{rowlen: N, data: a};
    let mut b = unsafe{ arena.alloc::<[f64;N]>() };

    // let mut a = Vec::with_capacity(N*N);
    // let mut b = Vec::with_capacity(N);
    // unsafe {
    //     // it's safe as long as uninitialized values aren't read
    //     a.set_len(N*N);
    //     b.set_len(N);
    // }
    // let mut a = A2DSlice{ rowlen: N, data: a.as_mut_slice()};
    
    for i in 0..N {
        for j in 0..N {
            a[I(i, j)] = rng.gen_range(0, RAND_MAX) as f64/10000.;
        }
        b[i] = 0.0f64;
    }

    for k in 0..N-1 {
        let m = (k+1..N).max_by_key(|&x| F64Ord(a[I(x, k)].abs()) ).unwrap();

        b[k] = (m+1) as f64;
        if m+2 != k {
            b[N-1] = -b[N-1];
        }

        a.swapi(I(m, k), I(k, k));
        let c = a[I(k, k)];

        if c == 0. { continue };

        for i in k+1..N {
            a[I(i,k)] = -a[I(i,k)]/c;
        }

        for j in k+1..N {
            let c = a[I(m, j)];
            a.swapi(I(m, j), I(k, j));
            if c != 0. {
                for i in k+1..N {
                    a[I(i,j)] += a[I(i,k)]*c;
                }
            }
        }
    }
    println!("{}", a[I(0,0)]);
    println!("{:?}", cmp(0.0f64.log2(), 0./0.0f64));
    let tend = precise_time_ns();
    let mut stderr = ::std::io::stderr();
    writeln!(stderr, "Elapsed {}ms", (tend-tstart)as f64/1_000_000.).unwrap();
}

/*
#include <cstdio>
#include <cstdlib>
#include <cmath>

int main() {
    const int N = 1200;
    double* A;
    double* B;
    double C;
    int m;
    
    A = new double[N*N];
    B = new double[N]; 
    
    //srand(time(0));
    //srand( (unsigned)time( NULL ) );
    
    for (int i  = 0; i < N; i++) {
        for (int j =  0; j < N; j++) {
            A[N*i+j] = (double)rand() / ((double)10000);
        }
        B[i] = 0.0;
    }

    // Сильно порезанный фрагмент кода:    
    for (int k = 0; k < N-1; k++ )  {
        m = (k + 1);
        
        for (int i = k + 1; i < N; i++ )   {
            if ( fabs(A[N*(i) + (k)]) > fabs(A[N*(m-1) + k]) ) m = (i + 1);
        }
        B[k] = m;
        
        if ( (m + 1) != k ) B[N-1] = -B[N-1];
        
        C = A[N*(m-1) + k];
        A[N*(m-1) + k] = A[N*k + k];
        A[N*k + k] = C;
        
        if (C == 0) continue;
        
        for (int i = k + 1; i < N; i++ ) A[N*i + k] = -A[N*i + k]/C;
        
        for (int j = k + 1; j < N; j++ )  {
            C = A[N*(m-1)+ j];
            A[N*(m-1) + j] = A[N*k + j];
            A[N*k + j] = C;
            if ( C != 0 )  {
                // дольше всего выполняется этот блок цикла:
                for (int i = k + 1; i < N; i++ )  {
                    A[N*i + j] = A[N*i + j] + A[N*i + k]*C;
                }
            }    
        }
    }
    
    delete [] A;
    delete [] B; 
    
    return 0;
}
*/
#[derive(Debug, Clone, Copy)]
struct I(usize, usize);

trait Array2D {
    type Item;
    fn data(&self) -> &[Self::Item];
    fn data_mut(&mut self) -> &mut [Self::Item];
    fn idx(&self, idx: I) -> usize;
}

struct A2DSlice<'a, T: 'a> {
    rowlen: usize,
    data: &'a mut [T],
}

impl<'a, T> Array2D for A2DSlice<'a, T> {
    type Item = T;
    fn data(&self) -> &[T] {
        self.data
    }
    fn data_mut(&mut self) -> &mut [T] {
        self.data
    }
    fn idx(&self, idx: I) -> usize {
        idx.0+idx.1*self.rowlen
    }
}

impl<'a, T> Index<I> for A2DSlice<'a, T> {
    type Output = T;
    fn index(&self, i: I) -> &Self::Output {
        let pos = self.idx(i);
        &self.data()[pos]
    }
}

impl<'a, T> IndexMut<I> for A2DSlice<'a, T> {
    fn index_mut(&mut self, i: I) -> &mut T {
        let pos = self.idx(i);
        &mut self.data_mut()[pos]
    }
}

trait SwapI: Array2D {
    fn swapi(&mut self, a: I, b: I);
}

impl<'a, T: Array2D> SwapI for T {
    #[inline]
    fn swapi(&mut self, a: I, b: I) {
        let a = self.idx(a);
        let b = self.idx(b);
        self.data_mut().swap(a, b);
    }
}

#[derive(Debug, Clone, Copy, PartialOrd)]
struct F64Ord(f64);

impl PartialEq for F64Ord {
    fn eq(&self, other: &Self) -> bool {
        if self.0.is_nan() && other.0.is_nan() {
            unsafe {
                let a: u64 = ::std::mem::transmute(self.0);
                let b: u64 = ::std::mem::transmute(other.0);
                a == b
            }
        } else {
            self.0 == other.0
        }
    }
}

impl Eq for F64Ord {}

impl Ord for F64Ord {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.partial_cmp(other) {
            Some(ord) => ord,
            None => {
                unsafe {
                    let a: u64 = ::std::mem::transmute(self.0);
                    let b: u64 = ::std::mem::transmute(other.0);
                    a.cmp(&b)
                }
            }
        }
    }
}

#[inline]
fn cmp(a: f64, b: f64) -> Ordering {
    match a.partial_cmp(&b) {
        Some(ord) => ord,
        None => {
            unsafe {
                let a: u64 = ::std::mem::transmute(a);
                let b: u64 = ::std::mem::transmute(b);
                a.cmp(&b)
            }
        }
    }
}

mod large_page_alloc {
    use advapi32::*;
    use winapi::*;
    use kernel32::*;
    use std::ptr::{null_mut};
    use std::cell::RefCell;

    fn error(err: &str) -> ! {
        panic!("{} Error 0x{:08x} ", err, unsafe{GetLastError()});
    }

    #[repr(C)]
    struct TKP {
        tkp: TOKEN_PRIVILEGES,
        laa: LUID_AND_ATTRIBUTES,
    }

    pub struct LPArena {
        inner: RefCell<InnerArena>,
    }

    #[derive(Debug, Clone, Copy)]
    struct InnerArena {
        base: usize,
        size: usize,
        ptr: usize,
    }


    pub fn init() -> LPArena {
        unsafe {
            let mut h_token = 0 as *mut _;
            if OpenProcessToken(GetCurrentProcess(), TOKEN_ADJUST_PRIVILEGES | TOKEN_QUERY, &mut h_token) == 0 {
                error("OpenProcessToken");
            }
            let mut tkp = ::std::mem::zeroed::<TKP>();
            assert_eq!(tkp.tkp.Privileges.len(), 0);
            let c = ::std::ffi::CString::new("SeLockMemoryPrivilege").unwrap();
            if LookupPrivilegeValueA(0 as *const _, c.as_ptr(), &mut tkp.laa.Luid) == 0 {
                error("LookupPrivilegeValueA");
            }
            tkp.tkp.PrivilegeCount = 1;
            tkp.laa.Attributes = SE_PRIVILEGE_ENABLED;
            if AdjustTokenPrivileges(h_token, FALSE, &mut tkp.tkp, 0, 0 as *mut _, 0 as *mut _) == 0 {
                error("AdjustTokenPrivileges");
            }
            let lpmin = GetLargePageMinimum();
            println!("Minimum large page size: {}", lpmin);
            if lpmin == 0 {
                error("GetLargePageMinimum");
            }
            let size = lpmin*40;
            let ptr = VirtualAlloc(null_mut(), size, MEM_RESERVE | MEM_COMMIT | MEM_LARGE_PAGES, PAGE_READWRITE);
            //let ptr = VirtualAlloc(null_mut(), size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
            if ptr == null_mut() {
                error("VirtualAlloc");
            }
            LPArena {
                inner: RefCell::new(InnerArena{
                    base: ptr as usize,
                    size: size as usize,
                    ptr: ptr as usize,
            })}
        }
    }

    impl LPArena {
        // returns reference to uninitialized value
        pub unsafe fn alloc<T: Sized>(&self) -> &mut T {
            let mut s = self.inner.borrow_mut();
            let size = ::std::mem::size_of::<T>();
            // panics on zero-sized types
            let next_addr = if s.ptr % size == 0 {s.ptr} else {((s.ptr/size)+1)*size}; 
            s.ptr = next_addr + size;
            if s.ptr >= s.base + s.size {
                panic!("Out of arena memory");
            }
            (next_addr as *mut _).as_mut().unwrap()
        }
    }

    impl Drop for LPArena {
        fn drop(&mut self) {
            let s = self.inner.borrow();
            unsafe{
                if VirtualFree(s.base as *mut _, 0, MEM_RELEASE) == 0 {
                    error("VirtualFree");
                }
            }
        }
    }
}

