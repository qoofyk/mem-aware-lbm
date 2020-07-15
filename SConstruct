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
step2_whole_Flags = argdict['step2_whole_Flags'].lower() == 'true'
step2_3parts_Flags = argdict['step2_3parts_Flags'].lower() == 'true'
step2_omp_Flags  = argdict['step2_omp_Flags'].lower() == 'true'
step2_unroll_Flags  = argdict['step2_unroll_Flags'].lower() == 'true'
step2_pyramid_Flags = argdict['step2_pyramid_Flags'].lower() == 'true'
if 'pillar_mem_Flags' in argdict:
    pillar_mem_Flags = argdict['pillar_mem_Flags'].lower() == 'true'
else: 
    pillar_mem_Flags = False
    
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

if step2_whole_Flags:
    flags.append('-DSTEP2_WHOLE')

if step2_3parts_Flags:
    flags.append('-DSTEP2_3PARTS')

if step2_omp_Flags:
    flags.append('-DSTEP2_OMP')

if step2_unroll_Flags:
    flags.append('-DSTEP2_UNROLL')

if step2_pyramid_Flags:
    flags.append('-DSTEP2_PYRAMID')

if pillar_mem_Flags:
    flags.append('-DPILLAR_MEM')

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
    if step2_whole_Flags:
        if step2_omp_Flags:
            if step2_unroll_Flags:
                if pillar_mem_Flags:
                    if '-DCUBE_MAP' in flags:
                        palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_whole_omp_unroll_pyramid_pillar_mem_cube_map',
                                          source  = sourceFiles )
                    if '-DPILLAR_SEQ_OMP' in flags:
                        palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_whole_omp_unroll_pyramid_pillar_mem_seq_omp',
                                          source  = sourceFiles )
                    else:
                        palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_whole_omp_unroll_pyramid_pillar_mem_pillar_map',
                                          source  = sourceFiles )
                else:
                    palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_whole_omp_unroll_pyramid',
                                          source  = sourceFiles )
            else:
                palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_whole_omp_pyramid',
                                      source  = sourceFiles )
        else: # sequential
            if step2_unroll_Flags:
                palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_whole_seq_unroll_pyramid',
                                      source  = sourceFiles )
            else:
                palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_whole_seq_pyramid',
                                      source  = sourceFiles )
    elif step2_3parts_Flags: # 3 parts
        if step2_omp_Flags:
            if step2_pyramid_Flags:
                palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_3parts_omp_pyramid',
                                      source  = sourceFiles )
            elif pillar_mem_Flags:
                palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_3parts_omp_line_pillar',
                                  source  = sourceFiles )
            else:
                palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_3parts_omp_line',
                                  source  = sourceFiles )
        else: # sequential
            if step2_pyramid_Flags:
                palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_3parts_seq_pyramid',
                                      source  = sourceFiles )
            else:
                palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi_step2_3parts_seq_line',
                                      source  = sourceFiles )
    else: #original
        palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb_mpi',
                                  source  = sourceFiles )
else:
    palabos_library = LibraryGen( target  = palabosRoot+'/lib/plb',
                                  source  = sourceFiles )

local_objects = env.Object(source = projectFiles)

all_objects = local_objects + palabos_library + libraries
print(all_objects)

env.Program(all_objects, LIBPATH=libraryPaths)
