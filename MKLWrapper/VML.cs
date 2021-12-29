
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
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "vdAbs")]
        internal static extern void _vmldAbs(int n, double [] a, double [] r);
        public static void Abs(double[] a, double[] r) => _vmldAbs(a.Length,a,r);
        /// <summary>
        /// r<-a + b
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "vdAdd")]
        internal static extern void _vmldAdd(int n, double [] a, double [] b, double [] r);
        public static void Add(double[] a, double[] b, double[] r) => _vmldAdd(a.Length, a, b, r);
        /// <summary>
        /// r <- a-b
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "vdSub")]
        internal  static extern void _vmldSub(int n, double [] a, double [] b,double [] r);
        public static void Sub(double[] a, double [] b, double[] r) => _vmldSub(a.Length, a, b, r);
        /// <summary>
        /// r <- 1.0/a
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "vdInv")]
        internal static extern void _vmldInv(int n, double [] a, double [] r);
        public static void Inv(double[] a, double[] r) => _vmldInv(a.Length, a, r);
        /// <summary>
        /// r<- sqrt(a)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl,EntryPoint = "_vdSqrt")]
        internal static extern void _vmldSqrt(int n, double []  a, double [] r);
        public static void Sqrt(double[] a, double[] r) => _vmldSqrt(a.Length, a, r);

        /// <summary>
        /// r<-1/sqrt(a)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "vdInvSqrt")]
        internal static extern void _vmldInvSqrt(int n, double [] a, double [] r);
        public static void invSqrt(double[] a, double[] r) => _vmldInvSqrt(a.Length, a, r);
        /// <summary>
        /// r<-a^2
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "vdSqr")]
        internal static extern void _vmldSqr(int n, double [] a, double [] r);
        public static void Sqr(double[] a, double[] r) => _vmldSqr(a.Length, a, r);

        /// <summary>
        /// r<-exp(a)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "vdExp")]
        internal static extern void _vmldExp(int n, double [] a, double [] r);
        public static void Exp(double[] a, double[] r) => _vmldExp(a.Length, a, r);
        /// <summary>
        /// r<-ln(x)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "vdLn")]
        internal static extern void _vmldLn(int n, double[] a, double[] r);
        public static void Ln(double[] a, double[] r) => _vmldLn(a.Length, a, r);
        /// <summary>
        /// r<-Cos(a)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl,EntryPoint ="vdCos")]
        internal  static extern void _vmldCos(int n, double [] a, double [] r);
        public static void Cos(double[] a, double[] r) => _vmldCos(a.Length, a, r);


        /// <summary>
        /// r<-Sin(a)
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl,EntryPoint ="vdSin")]
        internal  static extern void _vmldSin(int n, double [] a, double [] r);
        public static void Sin(double[] a, double[] r) => _vmldSin(a.Length, a, r);

        /// <summary>
        /// Error fucntion 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl,EntryPoint ="vdErf")]
        public static extern void _vmldErf(int n, double* a, double* r);
        /// <summary>
        /// NormCdf
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl,EntryPoint ="vdCdfNorm")]
        internal  static extern void _vmldCdfNorm(int n, double [] a, double [] r);
        public static void CdfNorm(double[] a, double[] r) => _vmldCdfNorm(a.Length, a, r);

        /// <summary>
        /// InvNormCdf
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="r"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl)]
        internal  static extern void _vmldCdfNormInv(int n, double [] a, double [] r);
        public static void CdfNormInv(double[] a, double[] r) => _vmldCdfNormInv(a.Length, a, r);

        /// <summary>
        /// y[i] = a[ia[i]]
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="ia"></param>
        /// <param name="y"></param>
        [DllImport(Controls.VML_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint ="vdPackV")]
        internal static extern void _vmldPackV(int n, double [] a, int [] ia, double [] y);
        public static void PackV(double[] a, int[] ia, double[] y) => _vmldPackV(a.Length, a, ia, y);
    }
}
