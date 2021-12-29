using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Linq;

namespace IntelMKL
{
    public unsafe class MKL
    {
        public enum CBLAS_ORDER
        {
            CblasRowMajor = 101,    /* row-major arrays */
            CblasColMajor = 102     /* column-major arrays */
        };

        public enum LAPACK_ORDER
        {
            LAPACKRowMajor = 101,    /* row-major arrays */
            LAPACKColMajor = 102     /* column-major arrays */
        };

        public enum CBLAS_TRANSPOSE
        {
            CblasNoTrans = 111,     /* trans='N' */
            CblasTrans = 112,       /* trans='T' */
            CblasConjTrans = 113    /* trans='C' */
        };

        public enum CBLAS_UPLO
        {
            CblasUpper = 121,       /* uplo ='U' */
            CblasLower = 122        /* uplo ='L' */
        };

        public enum CBLAS_DIAG
        {
            CblasNonUnit = 131,     /* diag ='N' */
            CblasUnit = 132         /* diag ='U' */
        };

        public enum CBLAS_SIDE
        {
            CblasLeft = 141,        /* side ='L' */
            CblasRight = 142        /* side ='R' */
        };

        public enum VMLAccuracy
        {
            LowAccuracy = 1, HighAccuracy = 2, EnhancedPerformanceAccuracy = 3
        }
        unsafe internal struct MKLVersion_
        {
            public int MajorVersion;
            public int MinorVersion;
            public int UpdateVersion;
            public sbyte* ProductStatus;
            public sbyte* Build;
            public sbyte* Processor;
            public sbyte* Platform;
        }


        public struct MKLVersion
        {
            public int MajorVersion;
            public int MinorVersion;
            public int UpdateVersion;
            public string ProductStatus;
            public string Build;
            public string Processor;
            public string Platform;
        }

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        static extern int mkl_get_version(ref MKLVersion_ version);
        public static MKLVersion GetMKLVersion()
        {
            var mklVer_ = new MKLVersion_();
            mkl_get_version(ref mklVer_);
            unsafe
            {
                return new MKLVersion
                {
                    MajorVersion = mklVer_.MajorVersion,
                    MinorVersion = mklVer_.MinorVersion,
                    UpdateVersion = mklVer_.UpdateVersion,
                    ProductStatus = new string(mklVer_.ProductStatus),
                    Build = new string(mklVer_.Build),
                    Processor = new string(mklVer_.Processor),
                    Platform = new string(mklVer_.Platform),
                };
            }
        }

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        internal static extern void mkl_free_buffers();
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void free_buffers() => mkl_free_buffers();

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        internal  static extern void mkl_set_num_threads(ref int n);

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        internal  static extern int mkl_get_max_threads();

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl,  ExactSpelling = true, SetLastError = false)]
        internal static extern void mkl_set_dynamic(ref int n);

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl,  ExactSpelling = true, SetLastError = false)]
        internal static extern int mkl_get_dynamic();

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        internal  static extern int mkl_set_threading_layer(ref int layer);
        /// <summary>
        /// Aligned allocators 
        /// </summary>
        /// <param name="alloc_size"></param>
        /// <param name="alignment"></param>
        /// <returns></returns>
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        internal static extern IntPtr mkl_malloc(ref IntPtr alloc_size, ref int alignment);

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        internal static extern void mkl_free(ref IntPtr a_ptr);

        public static IntPtr mklMalloc(int alloc_size, int alignment)
        {
            IntPtr p = new IntPtr(alloc_size);
            return mkl_malloc(ref p, ref alignment);
        }
        public static IntPtr mklMalloc(long alloc_size, int alignment)
        {
            IntPtr p = new IntPtr(alloc_size);
            return mkl_malloc(ref p, ref alignment);
        }

        public static void mklFree(ref IntPtr p) => mkl_free(ref p);
        /// <summary>
        /// BLAS-DDOT
        /// </summary>
        /// <param name="N"></param>
        /// <param name="x"></param>
        /// <param name="incX"></param>
        /// <param name="y"></param>
        /// <param name="incY"></param>
        /// <returns></returns>
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        internal  static extern double cblas_ddot(int N, double [] x, int incX, double [] y, int incY);
        /// <summary>
        /// L2 norm
        /// </summary>
        /// <param name="N"></param>
        /// <param name="x"></param>
        /// <param name="incX"></param>
        /// <param name="y"></param>
        /// <param name="incY"></param>
        /// <returns></returns>
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        public static extern double cblas_dnrm2(int N, double[] x, int incX);
        /// <summary>
        /// L1 norm
        /// </summary>
        /// <param name="N"></param>
        /// <param name="x"></param>
        /// <param name="incX"></param>
        /// <returns></returns>
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        public static extern double cblas_dasum(int N, double[] x, int incX);
        /// <summary>
        /// 
        /// </summary>
        /// <param name="matrix_layout"></param>
        /// <param name="trans"></param>
        /// <param name="m"></param>
        /// <param name="n"></param>
        /// <param name="nrhs"></param>
        /// <param name="a"></param>
        /// <param name="lda"></param>
        /// <param name="b"></param>
        /// <param name="ldb"></param>
        /// <returns></returns>
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        public static extern int LAPACKE_dgels(int matrix_layout, char trans, int m, int n, int nrhs, double[] a, int lda, double[] b, int ldb);

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        static extern double cblas_ddot(int N, double* X, int incX, double* Y, int incY);
        public static double dot(int N, double[] X, int iniX, int incX, double[] Y, int iniY, int incY)
        {
            fixed (double* xp = &X[iniX])
            fixed (double* yp = &Y[iniY])
            {
                return cblas_ddot(N, xp, incX, yp, incY);
            }
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double dot(double[] X, double[] Y) => dot(X.Length, X, 0, 1, Y, 0, 1);
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        static extern void cblas_dswap(int N, double* X, int incX, double* Y, int incY);
        public static void swap(int N, double[] X, int iniX, int incX, double[] Y, int iniY, int incY)
        {
            fixed (double* xp = &X[iniX])
            fixed (double* yp = &Y[iniY])
                cblas_dswap(N, xp, incX, yp, incY);
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void swap(double[] X, double[] Y) => swap(X.Length, X, 0, 1, Y, 0, 1);
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        static extern void cblas_dcopy(int N, double* X, int incX, double* Y, int incY);
        public static void copy(int N, double[] X, int iniX, int incX, double[] Y, int iniY, int incY)
        {
            fixed (double* xp = &X[iniX])
            fixed (double* yp = &Y[iniY])
                cblas_dcopy(N, xp, incX, yp, incY);
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void copy(double[] X, double[] Y) => copy(X.Length, X, 0, 1, Y, 0, 1);

        /// <summary>
        /// BLAS methods: daxpy
        /// </summary>
        /// <param name="N"></param>
        /// <param name="x"></param>
        /// <param name="incX"></param>
        /// <param name="y"></param>
        /// <param name="incY"></param>
        /// <returns></returns>
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        public static extern double cblas_daxpy(int N, double a, double* x, int incX, double* y, int incY);
        public static void axpy(int N, double a, double[] X, int iniX, int incX, double[] Y, int iniY, int incY)
        {
            fixed (double* xp = &X[iniX])
            fixed (double* yp = &Y[iniY])
                cblas_daxpy(N, a, xp, incX, yp, incY);
        }
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void axpy(double a, double[] X, double[] Y) => axpy(X.Length, a, X, 0, 1, Y, 0, 1);


        #region --optimization--
        /// <summary>
        ///  NonLinear LeastSquares results 
        /// </summary>
        public enum SolveResult
        {
            EXCEEDED_ITERATIONS = -1,
            DELTA_LESS_THAN_EPS0 = -2,
            F_NORM_2_LESS_THAN_EPS1 = -3,
            JACOBIAN_SINGULAR_LESS_THAN_ESP2 = -4,
            S_NORM_2_LESS_THAN_EPS3 = -5,
            F_J_S_LESS_THAN_EPS4 = -6,
            INVALID_OPTION = 1502,
            OUT_OF_MEMORY = 1503,
            NAN_PARAMETER = 10000,
        }
        //define the delegate 
        public unsafe delegate void SolveFn(int* m, int* n, double* x, double* fvec);
        //MKL entry points 
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        internal static extern int dtrnlsp_init(IntPtr* handle, int* n, int* m, double* x, double* eps, int* iter1, int* iter2, double* rs);
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        internal static extern int dtrnlsp_solve(IntPtr* handle, double* fvec, double* fjac, int* RCI_Request);
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        internal static extern int dtrnlsp_delete(IntPtr* handle);
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        internal static extern int dtrnlspbc_init(IntPtr* handle, int* n, int* m, double* x, double* lower, double* upper, double* eps, int* iter1, int* iter2, double* rs);
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        internal static extern int dtrnlspbc_solve(IntPtr* handle, double* fvec, double* fjac, int* RCI_Request);
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        internal static extern int dtrnlspbc_delete(IntPtr* handle);
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        internal static extern int djacobi_init(IntPtr* handle, int* n, int* m, double* x, double* fjac, double* eps);
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        internal static extern int djacobi_solve(IntPtr* handle, double* f1, double* f2, int* RCI_Request);
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        internal static extern int djacobi_delete(IntPtr* handle);
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true)]
        internal static extern int djacobi(SolveFn fcn, int* n, int* m, double* fjac, double* x, double* eps);


        const int ONE_ITERATION = 0;
        const int CALCULATE_FUNCTION = 1;
        const int CALCULATE_JACOBIAN = 2;
        const int SUCCESS = 1501;
        //default tolerance[6] array 
        public static readonly double[] DEFAULT_TOL = new[] { 1e-6, 1e-5, 1e-5, 1e-2, 1e-2, 1e-2 };
        public static bool IsNaN(double[] x)
        {
          return  x.Any(v => Double.IsNaN(v));
        }
        public static bool Isinf(double[] x)
        {
          return x.Any(v => Double.IsInfinity(v));
        }

        /// <summary>
        /// Nonlinear Least Squares Solver. Fx = F(x). minimizes Norm_2(F(x))
        /// </summary>
        /// <param name="F">objective function</param>
        /// <param name="x">input values, initial guess to start, solution on exit</param>
        /// <param name="Fx">function values, just zero to start, solution on exit</param>
        /// <param name="eps">precisions for stop-criteria, defaults to all 1e-9</param>
        /// <param name="iter1">maximum number of iterations, defaults of 1000</param>
        /// <param name="iter2">maximum number of iterations of calculation of trial-step, default of 100</param>
        /// <param name="rs">initial step bound (0.1 - 100.0 recommended), default of 0.0 which MKL defaults as 100.0</param>
        /// <param name="jeps">precision of the Jacobian matrix calculation</param>
        /// <returns>stop criterion</returns>
        public static SolveResult NonLinearLeastSquares(Action<double[], double[]> F, 
            double[] x, double[] Fx, double[] eps, 
            int iter1 = 100, int iter2 = 10, 
            double rs = 0.0, double jeps = 1e-7)
        {
            int n = x.Length;
            int m = Fx.Length;
            var jac = new double[n * m];
            var f1 = new double[m];
            var f2 = new double[m];
            fixed (double* xp = &x[0], epsp = &eps[0], Fxp = &Fx[0], jacp = &jac[0], f1p = &f1[0], f2p = &f2[0])
            {
                IntPtr handle;
                int request;
                int status = dtrnlsp_init(&handle, &n, &m, xp, epsp, &iter1, &iter2, &rs);

                while (status == SUCCESS && (status = dtrnlsp_solve(&handle, Fxp, jacp, &request)) == SUCCESS)
                {
                    if (request == CALCULATE_FUNCTION)
                    {
                        if (MKL.IsNaN(x)) return SolveResult.NAN_PARAMETER;
                        F(x, Fx);
                    }
                    else if (request == CALCULATE_JACOBIAN)
                    {
                        IntPtr jacHandle;
                        int jacRequest;
                        status = djacobi_init(&jacHandle, &n, &m, xp, jacp, &jeps);
                        if (status == SUCCESS)
                            while ((status = djacobi_solve(&jacHandle, f1p, f2p, &jacRequest)) == SUCCESS && jacRequest != 0)
                            {
                                if (MKL.IsNaN(x)) return SolveResult.NAN_PARAMETER;
                                F(x, jacRequest == 1 ? f1 : f2);
                            }
                        djacobi_delete(&jacHandle);
                    }
                    else if (request != ONE_ITERATION)
                    {
                        status = request;
                        break;
                    }
                }
                dtrnlsp_delete(&handle);
                return (SolveResult)status;
            }
        }

        /// <summary>
        /// Nonlinear Least Squares Solver. Fx = F(x). minimizes Norm_2(F(x))
        /// </summary>
        /// <param name="F">objective function</param>
        /// <param name="x">input values, initial guess to start, solution on exit</param>
        /// <param name="Fx">function values, just zero to start, solution on exit</param>
        /// <param name="eps">precisions for stop-criteria, defaults to all 1e-9</param>
        /// <param name="iter1">maximum number of iterations, defaults of 1000</param>
        /// <param name="iter2">maximum number of iterations of calculation of trial-step, default of 100</param>
        /// <param name="rs">initial step bound (0.1 - 100.0 recommended), default of 0.0 which MKL defaults as 100.0</param>
        /// <param name="jeps">precision of the Jacobian matrix calculation</param>
        /// <returns>stop criterion</returns>
        public static SolveResult NonLinearLeastSquares(SolveFn F, double[] x, double[] Fx, double[] eps,
                                int iter1 = 100, int iter2 = 10, double rs = 0.0, double jeps = 1e-6)
        {
            int n = x.Length;
            int m = Fx.Length;
            var jac = new double[n * m];
            fixed (double* xp = &x[0], epsp = &eps[0], Fxp = &Fx[0], jacp = &jac[0])
            {
                IntPtr handle;
                int request;
                var status = dtrnlsp_init(&handle, &n, &m, xp, epsp, &iter1, &iter2, &rs);
                while (status == SUCCESS && (status = dtrnlsp_solve(&handle, Fxp, jacp, &request)) == SUCCESS)
                {
                    if (request == CALCULATE_FUNCTION)
                    {
                        if (IsNaN(x)) return SolveResult.NAN_PARAMETER;
                        F(&m, &n, xp, Fxp);
                    }
                    else if (request == CALCULATE_JACOBIAN)
                    {
                        if (IsNaN(x)) return SolveResult.NAN_PARAMETER;
                        status = djacobi(F, &n, &m, jacp, xp, &jeps);
                    }
                    else if (request != ONE_ITERATION)
                    {
                        status = request;
                        break;
                    }
                }
                dtrnlsp_delete(&handle);
                return (SolveResult)status;
            }
        }


        /// <summary>
        /// Nonlinear Least Squares Solver. Fx = F(x). minimizes Norm_2(F(x))
        /// </summary>
        /// <param name="F">objective function</param>
        /// <param name="J">Jacobian function, J_ij = df_i / dx_j as a column major array</param>
        /// <param name="x">input values, initial guess to start, solution on exit</param>
        /// <param name="Fx">function values, just zero to start, solution on exit</param>
        /// <param name="eps">precisions for stop-criteria, defaults to all 1e-9</param>
        /// <param name="iter1">maximum number of iterations, defaults of 1000</param>
        /// <param name="iter2">maximum number of iterations of calculation of trial-step, default of 100</param>
        /// <param name="rs">initial step bound (0.1 - 100.0 recommended), default of 0.0 which MKL defaults as 100.0</param>
        /// <returns>stop criterion</returns>
        public static SolveResult NonLinearLeastSquares(SolveFn F, Action<double[], double[]> J,
                                                        double[] x, double[] Fx, double[] eps,
                                                        int iter1 = 100, int iter2 = 10, double rs = 0.0)
        {
            int n = x.Length;
            int m = Fx.Length;
            var jac = new double[n * m];
            fixed (double* xp = &x[0], epsp = &eps[0], Fxp = &Fx[0], jacp = &jac[0])
            {
                IntPtr handle;
                int request;
                var status = dtrnlsp_init(&handle, &n, &m, xp, epsp, &iter1, &iter2, &rs);
                if (status == SUCCESS)
                    while ((status = dtrnlsp_solve(&handle, Fxp, jacp, &request)) == SUCCESS)
                    {
                        if (request == CALCULATE_FUNCTION)
                        {
                            if (MKL.IsNaN(x)) return SolveResult.NAN_PARAMETER;
                            F(&m, &n , xp, Fxp);
                        }
                        else if (request == CALCULATE_JACOBIAN)
                        {
                            if (MKL.IsNaN(x)) return SolveResult.NAN_PARAMETER;
                            J(x, jac);
                        }
                        else if (request != ONE_ITERATION)
                        {
                            status = request;
                            break;
                        }
                    }
                dtrnlsp_delete(&handle);
                return (SolveResult)status;
            }
        }

        // <summary>
        /// Nonlinear Least Squares Solver. Fx = F(x). minimizes Norm_2(F(x))
        /// </summary>
        /// <param name="F">objective function, called in parallel</param>
        /// <param name="J">Jacobian function, J_ij = df_i / dx_j as a column major array</param>
        /// <param name="x">input values, initial guess to start, solution on exit</param>
        /// <param name="lower">x lower bound</param>
        /// <param name="upper">x upper bound</param>
        /// <param name="Fx">function values, just zero to start, solution on exit</param>
        /// <param name="eps">precisions for stop-criteria, defaults to all 1e-9</param>
        /// <param name="iter1">maximum number of iterations, defaults of 1000</param>
        /// <param name="iter2">maximum number of iterations of calculation of trial-step, default of 100</param>
        /// <param name="rs">initial step bound (0.1 - 100.0 recommended), default of 0.0 which MKL defaults as 100.0</param>
        /// <returns>stop criterion</returns>
        public static SolveResult NonLinearLeastSquares(SolveFn F, Action<double[], double[]> J,
                                                        double[] x, double[] lower, double[] upper, double[] Fx, double[] eps, 
                                                        int iter1 = 100, int iter2 = 10, double rs = 0.1)
        {
            int n = x.Length;
            int m = Fx.Length;
            var jac = new double[n * m];
            var f1 = new double[m];
            var f2 = new double[m];
            fixed (double* xp = &x[0], lowerp = &lower[0], upperp = &upper[0],Fxp = &Fx[0], epsp = &eps[0],
                jacp = &jac[0], f1p = &f1[0], f2p = &f2[0])
            {
                IntPtr handle;
                int request;
                var status = dtrnlspbc_init(&handle, &n, &m, xp, lowerp, upperp, epsp, &iter1, &iter2, &rs);
                if (status == SUCCESS)
                    while ((status = dtrnlspbc_solve(&handle, Fxp, jacp, &request)) == SUCCESS)
                    {
                        if (request == CALCULATE_FUNCTION)
                        {
                            if (IsNaN(x)) return SolveResult.NAN_PARAMETER;
                            F(&m, &n, xp, Fxp);
                        }
                        else if (request == CALCULATE_JACOBIAN)
                        {
                            if (IsNaN(x)) return SolveResult.NAN_PARAMETER;
                            J(x, jac);
                        }
                        else if (request != ONE_ITERATION)
                        {
                            status = request;
                            break;
                        }
                    }
                dtrnlspbc_delete(&handle);
                return (SolveResult)status;
            }
        }

        /// <summary>
        /// Nonlinear Least Squares Solver. Fx = F(x). minimizes Norm_2(F(x))
        /// </summary>
        /// <param name="F">objective function, not called in parallel</param>
        /// <param name="x">input values, initial guess to start, solution on exit</param>
        /// <param name="lower">x lower bound</param>
        /// <param name="upper">x upper bound</param>
        /// <param name="Fx">function values, just zero to start, solution on exit</param>
        /// <param name="eps">precisions for stop-criteria, defaults to all 1e-9</param>
        /// <param name="iter1">maximum number of iterations, defaults of 1000</param>
        /// <param name="iter2">maximum number of iterations of calculation of trial-step, default of 100</param>
        /// <param name="rs">initial step bound (0.1 - 100.0 recommended), default of 0.0 which MKL defaults as 100.0</param>
        /// <param name="jeps">precision of the Jacobian matrix calculation</param>
        /// <returns>stop criterion</returns>
        public static SolveResult NonLinearLeastSquaresBounded(Action<double[], double[]> F,
                double[] x, double[] lower, double[] upper, double[] Fx, double[] eps, 
                int iter1 = 100, int iter2 = 10, double rs = 0.1, double jeps = 1e-7)
        {
            //1. dtrnlsp_init(ref handle, ref n, ref m, x, eps, ref iter1, ref iter2, ref rs)
            //2.dtrnlsp_check(ref handle, ref n, ref m, fjac, fvec, eps, info)
            //3.LOOP dtrnlsp_solve(ref handle, fvec, fjac, ref RCI_Request)
            //4.dtrnlsp_get(ref handle, ref iter, ref st_cr, ref r1, ref r2)
           //5.dtrnlsp_delete(ref handle)

            int n = x.Length;
            int m = Fx.Length;
            var jac = new double[n * m];
            var f1 = new double[m];
            var f2 = new double[m];
            fixed (double* xp = &x[0], lowerp = &lower[0], upperp = &upper[0], epsp = &eps[0], Fxp = &Fx[0],
                jacp = &jac[0], f1p = &f1[0], f2p = &f2[0])
            {
                IntPtr handle;
                int RCI_request;
                int status = dtrnlspbc_init(&handle, &n, &m, xp, lowerp, upperp, epsp, &iter1, &iter2, &rs);
                while (status == SUCCESS && (status = dtrnlspbc_solve(&handle, Fxp, jacp, &RCI_request)) == SUCCESS)
                {
                    if (RCI_request == CALCULATE_FUNCTION)
                    {
                        if (MKL.IsNaN(x)) return SolveResult.NAN_PARAMETER;
                        F(x, Fx);
                    }
                    else if (RCI_request == CALCULATE_JACOBIAN)
                    {
                        IntPtr jacHandle;
                        int jacRequest;
                        status = djacobi_init(&jacHandle, &n, &m, xp, jacp, &jeps);
                        if (status == SUCCESS)
                            while ((status = djacobi_solve(&jacHandle, f1p, f2p, &jacRequest)) == SUCCESS && jacRequest != 0)
                            {
                                if (MKL.IsNaN(x)) return SolveResult.NAN_PARAMETER;
                                F(x, jacRequest == 1 ? f1 : f2);
                            }
                        djacobi_delete(&jacHandle);
                    }
                    else if (RCI_request != ONE_ITERATION)
                    {
                        status = RCI_request;
                        break;
                    }
                }
                dtrnlspbc_delete(&handle);
                return (SolveResult)status;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="F"></param>
        /// <param name="J"></param>
        /// <param name="x"></param>
        /// <param name="Fx"></param>
        /// <param name="eps"></param>
        /// <param name="iter1"></param>
        /// <param name="iter2"></param>
        /// <param name="rs"></param>
        /// <returns></returns>
        public static SolveResult NonLinearLeastSquaresBounded(Action<double[], double[]> F, Action<double[], double[]> J, double[] x, double[] Fx, double[] eps, int iter1 = 1000, int iter2 = 100, double rs = 0.0f)
        {
            int n = x.Length;
            int m = Fx.Length;
            var jac = new double[n * m];
            fixed (double* xp = &x[0], epsp = &eps[0], Fxp = &Fx[0], jacp = &jac[0])
            {
                IntPtr handle;
                int request;
                var status = dtrnlsp_init(&handle, &n, &m, xp, epsp, &iter1, &iter2, &rs);
                if (status == SUCCESS)
                    while ((status = dtrnlsp_solve(&handle, Fxp, jacp, &request)) == SUCCESS)
                    {
                        if (request == CALCULATE_FUNCTION)
                        {
                            if (IsNaN(x)) return SolveResult.NAN_PARAMETER;
                            F(x, Fx);
                        }
                        else if (request == CALCULATE_JACOBIAN)
                        {
                            if (IsNaN(x)) return SolveResult.NAN_PARAMETER;
                            J(x, jac);
                        }
                        else if (request != ONE_ITERATION)
                        {
                            status = request;
                            break;
                        }
                    }
                dtrnlsp_delete(&handle);
                return (SolveResult)status;
            }
        }

        #endregion
    }
}
