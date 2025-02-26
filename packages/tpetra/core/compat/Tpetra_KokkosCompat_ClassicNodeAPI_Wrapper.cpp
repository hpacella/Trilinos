#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include "Teuchos_ParameterList.hpp"

namespace Tpetra {
namespace KokkosCompat {

#ifdef KOKKOS_ENABLE_THREADS
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Threads>::name () {
      return "Threads/Wrapper";
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::OpenMP>::name () {
      return "OpenMP/Wrapper";
    }
#endif

#ifdef KOKKOS_ENABLE_SERIAL
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Serial>::name () {
      return "Serial/Wrapper";
    }
#endif // KOKKOS_ENABLE_SERIAL

#ifdef KOKKOS_ENABLE_CUDA
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Cuda>::name() {
      return std::string("Cuda/Wrapper");
    }
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_HIP
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Experimental::HIP, Kokkos::Experimental::HIPSpace>::name() {
      return std::string("HIP/Wrapper");
    }
#endif // KOKKOS_ENABLE_HIP

#ifdef KOKKOS_ENABLE_SYCL
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace>::name() {
      return std::string("SYCL/Wrapper");
    }
    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLSharedUSMSpace>::name() {
      return std::string("SYCL/Wrapper");
    }
#endif // KOKKOS_ENABLE_SYCL


} // namespace KokkosCompat
} // namespace Tpetra



