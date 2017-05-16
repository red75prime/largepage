extern crate rand;
extern crate advapi32;
extern crate kernel32;
extern crate winapi;
use rand::{StdRng, Rng};
use std::ops::{Index, IndexMut};
use std::cmp::Ordering;

mod large_page_alloc {
    use advapi32::*;
    use winapi::*;
    use kernel32::*;
    use std::ptr::{null_mut};

    fn error(err: &str) -> ! {
        panic!("{} Error 0x{:08x} ", err, unsafe{GetLastError()});
    }

    #[repr(C)]
    struct TKP {
        tkp: TOKEN_PRIVILEGES,
        laa: LUID_AND_ATTRIBUTES,
    }

    pub struct LPArena {
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
            if ptr == null_mut() {
                error("VirtualAlloc");
            }
            LPArena {
                base: ptr as usize,
                size: size as usize,
                ptr: ptr as usize,
            }
        }
    }

    impl LPArena {
        pub fn alloc<T: Sized>(&mut self) -> &'static mut T {
            let size = ::std::mem::size_of::<T>();
            // panics on zero-sized types
            let next_addr = if self.ptr%size == 0 {self.ptr} else {((self.ptr/size)+1)*size}; 
            self.ptr = next_addr + size;
            if self.ptr >= self.base + self.size {
                panic!("Out of arena memory");
            }
            unsafe {
                (next_addr as *mut _).as_mut().unwrap()
            }
        }
    }
}

const RAND_MAX: u32 = 32768;

const N: usize = 1200;

#[derive(Debug, Clone, Copy)]
struct I(usize, usize);

impl I {
    #[inline]
    fn idx(self) -> usize {
        self.0*N+self.1
    }
}

impl Index<I> for [f64] {
    type Output = f64;
    fn index(&self, i: I) -> &f64 {
        &self[i.idx()]
    }
}

impl IndexMut<I> for [f64] {
    //type Output = f64;
    fn index_mut(&mut self, i: I) -> &mut f64 {
        &mut self[i.idx()]
    }
}

impl Index<I> for Vec<f64> {
    type Output = f64;
    fn index(&self, i: I) -> &f64 {
        &self[i.idx()]
    }
}

impl IndexMut<I> for Vec<f64> {
    //type Output = f64;
    fn index_mut(&mut self, i: I) -> &mut f64 {
        &mut self[i.idx()]
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

fn main() {
    let mut arena = large_page_alloc::init();

    let mut rng = StdRng::new().expect("Cannot create rng");


//    let mut a = arena.alloc::<[f64;N*N]>();
//    let mut b = arena.alloc::<[f64;N]>();
    let mut a = Vec::with_capacity(N*N);
    let mut b = Vec::with_capacity(N);
    unsafe {
        // it's safe as long as uninitialized values aren't read
        a.set_len(N*N);
        b.set_len(N);
    }
    
    for i in 0..N {
        for j in 0..N {
            a[I(i, j)] = rng.gen_range(0, RAND_MAX) as f64/10000.;
        }
        b[i] = 0.0f64;
    }

    for k in 0..N-1 {
        let m = (k+1..N).max_by(|&x, &y| cmp(a[I(x, k)].abs(), a[I(y, k)].abs())).unwrap()+1;

        b[k] = m as f64;
        if m+1 != k {
            b[N-1] = -b[N-1];
        }

        a.swap(I(m-1, k).idx(), I(k, k).idx());
        let c = a[I(k, k)];

        if c == 0. { continue };

        for i in k+1..N {
            a[I(i,k)] = -a[I(i,k)]/c;
        }

        for j in k+1..N {
            let c = a[I(m-1, j)];
            a.swap(I(m-1, j).idx(), I(k, j).idx());
            if c != 0. {
                for i in k+1..N {
                    a[I(i,j)] += a[I(i,k)]*c;
                }
            }
        }
    }
    println!("{}", a[0]);
    println!("{:?}", cmp(0.0f64.log2(), 0./0.0f64))
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
