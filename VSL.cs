

using System.Security;
using System.Runtime.InteropServices;
using System;

namespace IntelMKL
{


    [SuppressUnmanagedCodeSecurity]
    public static unsafe class VSL
    {
        private const int VSL_BRNG_SHIFT = 20;
        private const int VSL_BRNG_INC = 1 << VSL_BRNG_SHIFT;
        public enum BRNG
        {
            VSL_BRNG_MCG31 = VSL_BRNG_INC,
            VSL_BRNG_R250 = (VSL_BRNG_MCG31 + VSL_BRNG_INC),
            VSL_BRNG_MRG32K3A = (VSL_BRNG_R250 + VSL_BRNG_INC),
            VSL_BRNG_MCG59    =  (VSL_BRNG_MRG32K3A +VSL_BRNG_INC),
            VSL_BRNG_WH      =   (VSL_BRNG_MCG59    +VSL_BRNG_INC),
             VSL_BRNG_SOBOL     = (VSL_BRNG_WH       +VSL_BRNG_INC),
             VSL_BRNG_NIEDERR    =(VSL_BRNG_SOBOL    +VSL_BRNG_INC),
             VSL_BRNG_MT19937    =(VSL_BRNG_NIEDERR  +VSL_BRNG_INC),
             VSL_BRNG_MT2203     =(VSL_BRNG_MT19937  +VSL_BRNG_INC),
             VSL_BRNG_IABSTRACT  =(VSL_BRNG_MT2203   +VSL_BRNG_INC),
             VSL_BRNG_DABSTRACT  =(VSL_BRNG_IABSTRACT+VSL_BRNG_INC),
             VSL_BRNG_SABSTRACT = (VSL_BRNG_DABSTRACT+VSL_BRNG_INC),
             VSL_BRNG_SFMT19937 = (VSL_BRNG_SABSTRACT+VSL_BRNG_INC),
             VSL_BRNG_NONDETERM = (VSL_BRNG_SFMT19937+VSL_BRNG_INC)

        }
        public enum InitMethod
        {
             VSL_INIT_METHOD_STANDARD = 0,
             VSL_INIT_METHOD_LEAPFROG  =1,
             VSL_INIT_METHOD_SKIPAHEAD =2

        }
        public enum GaussinMethod
        {
             VSL_RNG_METHOD_GAUSSIAN_BOXMULLER   =0, /* vsl{d,s}RngGaussian */
             VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2  =1, /* vsl{d,s}RngGaussian */
             VSL_RNG_METHOD_GAUSSIAN_ICDF       = 2 /* vsl{d,s}RngGaussian */
        }
        public enum ReturnCode
        {
            VSL_STATUS_OK =0,
            VSL_ERROR_FEATURE_NOT_IMPLEMENTED = -1,
            VSL_ERROR_UNKNOWN        =          -2
        }


        public class IntelMKLRandomNumberStream : IDisposable
        {
            [DllImport(Controls.VSL_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "vslNewStream", 
                ExactSpelling = true, SetLastError = false)]
            internal static extern int vslNewStream(ref IntPtr stream, int brng, int seed);

            [DllImport(Controls.VSL_DLL, CallingConvention = CallingConvention.Cdecl,
                ExactSpelling = true, SetLastError = false)]
            internal static extern int vslNewStreamEx(ref IntPtr stream, int brng, int nparams, uint[] parameters);

            [DllImport(Controls.VSL_DLL, CallingConvention = CallingConvention.Cdecl,
                ExactSpelling = true, SetLastError = false)]
            internal static extern int vslDeleteStream(ref IntPtr stream);

            public readonly IntPtr RandomStream;

            public IntelMKLRandomNumberStream(VSL.BRNG type, int seed)
            {
                RandomStream = new IntPtr();
                int status = -1;

                status = vslNewStream(ref RandomStream, (int)type, seed);
                //status = vslNewStream(ref randomStream, BRNGMapper[(int)BRNG.MRG32K3A], seed);

                if (status != 0)
                    throw new Exception("Random number generation stream failed to initialise.");
            }

            public void Dispose()
            {
                DeleteStream();
            }

            public static IntPtr[] CreateMultithreadStreams(int seed, int nStreams)
            {
                IntPtr[] multithreadStreams = new IntPtr[nStreams];
                for (int i = 0; i < nStreams; ++i)
                {
                    int status = -1;
                    status = vslNewStream(ref multithreadStreams[i], (int)BRNG.VSL_BRNG_MT19937, seed);
                    if (status != 0) throw new Exception("Random number generation stream failed to initialise.");
                }

                return multithreadStreams;
            }

            public void DeleteStream()
            {
                if (RandomStream != IntPtr.Zero)
                {
                    _DeleteStream(RandomStream);
                }
            }

            private void _DeleteStream(IntPtr stream)
            {
                int status = vslDeleteStream(ref stream);
                if (status != 0) throw new Exception("failed to delete stream.");
            }
            public static void DeleteStreams(IntPtr [] streams)
            {
                for(int i=0;i<streams.Length;++i)
                {
                    int status = vslDeleteStream(ref streams[i]);
                    if (status != 0) throw new Exception($"failed to delete stream# {i}.");
                    streams[i] = IntPtr.Zero;
                }
            }
        }

        #region random_number_generators
        /// <summary>
        /// 
        /// </summary>
        /// <param name="method"></param>
        /// <param name="stream"></param>
        /// <param name="length"></param>
        /// <param name="vector"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        [DllImport(Controls.VSL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern int vdRngUniform(int method, IntPtr stream, int length, [In, Out] double[] vector, double a, double b);
        /// <summary>
        /// 
        /// </summary>
        /// <param name="method"></param>
        /// <param name="stream"></param>
        /// <param name="length"></param>
        /// <param name="single"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        [DllImport(Controls.VSL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern int vdRngUniform(int method, IntPtr stream, int length, [In, Out] ref double single, double a, double b);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="method"></param>
        /// <param name="stream"></param>
        /// <param name="length"></param>
        /// <param name="vector"></param>
        /// <param name="mean"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        [DllImport(Controls.VSL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern int vdRngGaussian(int method, IntPtr stream, int length, double* vector, double mean, double sigma);
        /// <summary>
        /// 
        /// </summary>
        /// <param name="method"></param>
        /// <param name="stream"></param>
        /// <param name="length"></param>
        /// <param name="vector"></param>
        /// <param name="mean"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        [DllImport(Controls.VSL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern int vdRngGaussian(int method, IntPtr stream, int length, [In, Out] double[] vector, double mean, double sigma);
        /// <summary>
        /// 
        /// </summary>
        /// <param name="method"></param>
        /// <param name="stream"></param>
        /// <param name="length"></param>
        /// <param name="vector"></param>
        /// <param name="displacement"></param>
        /// <param name="rate"></param>
        /// <returns></returns>
        [DllImport(Controls.VSL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern int vdRngExponential(int method, IntPtr stream, int length, [In, Out] double[] vector, double displacement, double rate);
        #endregion
        #region threading_controls
        /// <summary>
        /// 
        /// </summary>
        /// <param name="nThreads"></param>
        /// <returns></returns>
        [DllImport(Controls.VSL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern int MKL_Set_Num_Threads(int nThreads);
        /// <summary>
        /// 
        /// </summary>
        /// <param name="threading"></param>
        /// <returns></returns>
        [DllImport(Controls.VSL_DLL, CallingConvention = CallingConvention.Cdecl,
            ExactSpelling = true, SetLastError = false)]
        public static extern int MKL_Set_Threading_Layer(int threading);
        #endregion
    }
}
