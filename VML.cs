
using System.Security;
using System.Runtime.InteropServices;

namespace IntelMKL
{
    [SuppressUnmanagedCodeSecurity]
    public static unsafe  class VML
    {
        /// <summary>
        /// r<- |a|
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_vmldAbs")]
        public static extern void _vmldAbs(int n, double [] a, double [] r);
        /// <summary>
        /// r< a+b
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_vmldAdd")]
        public static extern void _vmldAdd(int n, double [] a, double [] b, double [] r);
        /// <summary>
        /// r <- a-b
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_vmldSub")]
        public static extern void _vmldSub(int n, double [] a, double [] b,double [] r);
        /// <summary>
        /// r <- 1.0/a
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_vmldInv")]
        public static extern void _vmldInv(int n, double [] a, double [] r);
        /// <summary>
        /// r<- sqrt(a)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl,EntryPoint = "_vmldSqrt")]
        public static extern void _vmldSqrt(int n, double []  a, double [] r);

        /// <summary>
        /// r<-1/sqrt(a)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_vmldInvSqrt")]
        public static extern void _vmldInvSqrt(int n, double [] a, double [] r);
        /// <summary>
        /// r<-a^2
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_vmldSqr")]
        public static extern void _vmldSqr(int n, double [] a, double [] r);

        /// <summary>
        /// r<-exp(a)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_vmldExp")]
        public static extern void _vmldExp(int n, double [] a, double [] r);
        /// <summary>
        /// r<-ln(x)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_vmldLn")]
        public static extern void _vmldLn(int n, double[] a, double[] r);
        /// <summary>
        /// r<-Cos(a)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void _vmldCos(int n, double [] a, double [] r);
        /// <summary>
        /// r<-Sin(a)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void _vmldSin(int n, double [] a, double [] r);

        /// <summary>
        /// Error fucntion 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void _vmldErf(int n, double* a, double* r);
        /// <summary>
        /// NormCdf
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void _vmldCdfNorm(int n, double [] a, double [] r);

        /// <summary>
        /// InvNormCdf
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void _vmldCdfNormInv(int n, double [] a, double [] r);

        /// <summary>
        /// y[i] = a[ia[i]]
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="ia"></param>
        /// <param name="y"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void _vmldPackV(int n, double [] a, int [] ia, double [] y);
    }
}
