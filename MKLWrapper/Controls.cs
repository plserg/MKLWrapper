using System;
using System.Collections.Generic;
using System.Text;

namespace IntelMKL
{
    public static class Controls
    {
        /// <summary>
        ///  name of the oneAPI/MKL runtime library
        /// </summary>
        public const string MKL_DLL = "mkl_rt.2.dll";
        /// <summary>
        /// oneAPI/SIMD vector library 
        /// </summary>
        public const string VML_DLL = "mkl_avx2.2.dll";
        /// <summary>
        /// name of the stats MKL 
        /// </summary>
        public const string VSL_DLL = "mkl_rt.2.dll";
    }
}
