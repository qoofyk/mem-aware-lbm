###########################################################
# Configuration file for the compilation of Palabos code,
# using the SConstruct library.
# IT IS NOT RECOMMENDED TO MODIFY THIS FILE.
# Compilation should be personalized by adjusting the
# Makefile in the directory of the main source files.
# See Palabos examples for sample Makefiles.
###########################################################

import os
import sys
import glob

argdict = dict(ARGLIST)

# Read input parameters
palabosRoot   = argdict['palabosRoot']
projectFiles  = Split(argdict['projectFiles'])
optimize      = argdict['optimize'].lower() == 'true'
debug         = argdict['debug'].lower() == 'true'
profile       = argdict['profile'].lower() == 'true'
MPIparallel   = argdict['MPIparallel'].lower() == 'true'
SMPparallel   = argdict['SMPparallel'].lower() == 'true'
usePOSIX      = argdict['usePOSIX'].lower() == 'true'
serialCXX     = argdict['serialCXX']
parallelCXX   = argdict['parallelCXX']
compileFlags  = Split(argdict['compileFlags'])
linkFlags     = Split(argdict['linkFlags'])
optimFlags    = Split(argdict['optimFlags'])
debugFlags    = Split(argdict['debugFlags'])
profileFlags  = Split(argdict['profileFlags'])
libraryPaths  = Split(argdict['libraryPaths'])
includePaths  = Split(argdict['includePaths'])
libraries     = Split(argdict['libraries'])
step2_omp_Flags  = argdict['step2_omp_Flags'].lower() == 'true'
step2_seq_Flags  = argdict['step2_seq_Flags'].lower() == 'true'

# Read the optional input parameters
try:
    dynamicLibrary = argdict['dynamicLibrary'].lower() == 'true'
except:
    dynamicLibrary = False

try:
    srcPaths = Split(argdict['srcPaths'])
except:
    srcPaths = []

flags = compileFlags
allPaths = [palabosRoot+'/src'] + [palabosRoot+'/externalLibraries'] + includePaths

if optimize:
    flags.append(optimFlags)

if debug:
    flags.append(debugFlags)
    flags.append('-DPLB_DEBUG')

if profile:
    flags.append(profileFlags)
    linkFlags.append(profileFlags)

if MPIparallel:
    compiler = parallelCXX
    flags.append('-DPLB_MPI_PARALLEL')
else:
    compiler = serialCXX

if SMPparallel:
    flags.append('-DPLB_SMP_PARALLEL')

if usePOSIX:
    flags.append('-DPLB_USE_POSIX')

if step2_omp_Flags:
    flags.append('-DSTEP2_OMP')

env = Environment ( ENV       = os.environ,
                    CXX       = compiler,
                    CXXFLAGS  = flags,
                    LINKFLAGS = linkFlags,
                    CPPPATH   = allPaths
                  )

if dynamicLibrary:
    LibraryGen = env.SharedLibrary
else:
    LibraryGen = env.Library


sourceFiles = []
for srcDir in glob.glob(palabosRoot+'/src/*'):
    sourceFiles.extend(glob.glob(srcDir+'/*.cpp'))

for srcDir in srcPaths:
    sourceFiles.extend(glob.glob(srcDir+'/*.cpp'))

sourceFiles.extend(glob.glob(palabosRoot+'/externalLibraries/tinyxml/*.cpp'));

if MPIparallel:
    if step2_omp_Flags:
        palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_omp',
                                  source  = sourceFiles )
    elif step2_seq_Flags:
        palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_seq',
                                  source  = sourceFiles )
    else:
        palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi',
                                  source  = sourceFiles )
else:
    palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb',
                                  source  = sourceFiles )

local_objects = env.Object(source = projectFiles)

all_objects = local_objects + palabos_library

env.Program(all_objects, LIBS=libraries, LIBPATH=libraryPaths)
