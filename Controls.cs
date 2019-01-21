using System;
using System.Collections.Generic;
using System.Text;

namespace IntelMKL
{
    public static class Controls
    {
        /// <summary>
        ///  name of the MKL runtime library
        /// </summary>
        public const string MKL_DLL = "mkl_rt.dll";
        /// <summary>
        /// SIMD vector library 
        /// </summary>
        public const string VML_DLL = "mkl_vml_avx2.dll";
        /// <summary>
        /// name of the stats MKL 
        /// </summary>
        public const string VSL_DLL = "mkl_rt.dll";
    }
}
