#ifndef BLOCK_LATTICE_3D_PILLAR_MEM_H
#define BLOCK_LATTICE_3D_PILLAR_MEM_H
namespace plb {

// Add by Yuankun
#define YK_PANEL_LEN 32 // 2^5=5
#define YK_LOG_PANEL_LEN 5 // 2^5=5
#define YK_LOG_NY 9// 2^9 512
// #define YK_PANEL_LEN_MASK (~0 << YK_LOG_PANEL_LEN) //0xffffffe0
// End add by Yuankun

extern plint ykTile;

extern inline plint cube_mem_map_iX (plint iX, plint iY, plint iZ, plint nx, plint ny, plint nz);
extern inline plint pillar_mem_map_iX (plint iX, plint iY, plint iZ, plint nx, plint ny, plint nz);
}

#endif  // BLOCK_LATTICE_3D_PILLAR_MEM_H