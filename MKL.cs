using System;
using System.Collections.Generic;
using System.Text;
using System.Runtime.InteropServices;

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
        /// <summary>
        /// 
        /// </summary>
        /// <param name="N"></param>
        /// <param name="x"></param>
        /// <param name="incX"></param>
        /// <param name="y"></param>
        /// <param name="incY"></param>
        /// <returns></returns>
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        public static extern double cblas_daxpy(int N, double[] x, int incX, double[] y, int incY);
        /// <summary>
        /// 
        /// </summary>
        /// <param name="N"></param>
        /// <param name="x"></param>
        /// <param name="incX"></param>
        /// <param name="y"></param>
        /// <param name="incY"></param>
        /// <returns></returns>
        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl, ExactSpelling = true, SetLastError = false)]
        public static extern double cblas_ddot(int N, double [] x, int incX, double [] y, int incY);
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

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern void mkl_set_num_threads(ref int n);

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern int mkl_get_max_threads();

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern void mkl_set_dynamic(ref int n);

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern int mkl_get_dynamic();

        [DllImport(Controls.MKL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern int mkl_set_threading_layer(ref int layer);

    }
}
