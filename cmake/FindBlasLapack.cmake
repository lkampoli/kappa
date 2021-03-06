#
# this module look for OpenBlas (Blas, Lapack) support
# it will define the following values
#
# OPENBLAS_LIBRARY   = the library to link against (openblas etc)
# KAPPA_HAS_OPENBLAS = set to true after finding the library
#

# BLAS ############################

# FIND_LIBRARY(BLAS_LIBRARY blas ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/lib )
#
# IF ( BLAS_LIBRARY )
#  SET ( KAPPA_HAVE_BLAS 1 CACHE BOOL "Found BLAS library")
# ELSE()
#  SET ( KAPPA_HAVE_BLAS 0 )
# ENDIF()
#
# MESSAGE ( "KAPPA_HAVE_BLAS: [${KAPPA_HAVE_BLAS}]" )
# IF(KAPPA_HAVE_BLAS)
#  MESSAGE ( "  BLAS_LIBRARY:   [${BLAS_LIBRARY}]" )
# ENDIF(KAPPA_HAVE_BLAS)

# LAPACK #########################

# FIND_LIBRARY(LAPACK_LIBRARY lapack ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/lib)
#
# IF ( LAPACK_LIBRARY )
#  SET ( KAPPA_HAVE_LAPACK 1 CACHE BOOL "Found LAPACK library")
# ELSE()
#  SET ( KAPPA_HAVE_LAPACK 0 )
# ENDIF()
#
# MESSAGE ( "KAPPA_HAVE_LAPACK: [${KAPPA_HAVE_LAPACK}]" )
# IF(KAPPA_HAVE_LAPACK)
#  MESSAGE ( "  LAPACK_LIBRARY:   [${LAPACK_LIBRARY}]" )
# ENDIF(KAPPA_HAVE_LAPACK)

# OPENBLAS #########################

 FIND_LIBRARY(OPENBLAS_LIBRARY openblas ${CMAKE_CURRENT_SOURCE_DIR}/kappa_deps/x86_64/lib)

 IF ( OPENBLAS_LIBRARY )
  SET ( KAPPA_HAS_OPENBLAS 1 CACHE BOOL "Found OPENBLAS library")
 ELSE()
  SET ( KAPPA_HAS_OPENBLAS 0 )
 ENDIF()

 MESSAGE ( " KAPPA_HAS_OPENBLAS: [${KAPPA_HAS_OPENBLAS}]" )
 IF(KAPPA_HAS_OPENBLAS)
  MESSAGE ( " OPENBLAS_LIBRARY:  [${OPENBLAS_LIBRARY}]"   )
 ENDIF(KAPPA_HAS_OPENBLAS)

# BOTH ###########################

# IF ( CF_HAVE_BLAS AND CF_HAVE_LAPACK )
#  SET ( CF_BLASLAPACK_LIBRARIES   "${LAPACK_LIBRARY} ${BLAS_LIBRARY}" CACHE STRING "BLAS and LAPACK libraries")
#  SET ( CF_HAVE_BLASLAPACK ON CACHE BOOL "Found BLAS and LAPACK libraries")
# ENDIF()

 MARK_AS_ADVANCED ( LAPACK_LIBRARY BLAS_LIBRARY OPENBLAS_LIBRARY )

#################################

