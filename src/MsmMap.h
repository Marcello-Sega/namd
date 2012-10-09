/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef MSMMAP_H
#define MSMMAP_H

#define MSM_MIGRATION
#undef MSM_MIGRATION

#define MSM_MAX_BLOCK_SIZE 8
#define MSM_MAX_BLOCK_VOLUME \
  (MSM_MAX_BLOCK_SIZE * MSM_MAX_BLOCK_SIZE * MSM_MAX_BLOCK_SIZE)

#define DEBUG_MSM
#undef DEBUG_MSM

#define DEBUG_MSM_VERBOSE
#undef DEBUG_MSM_VERBOSE

#define DEBUG_MSM_GRID
#undef DEBUG_MSM_GRID

// assert macro
#undef ASSERT
#ifdef DEBUG_MSM
#define ASSERT(expr) \
  do { \
    if ( !(expr) ) { \
      char msg[100]; \
      snprintf(msg, sizeof(msg), "ASSERT: \"%s\" " \
          "(%s, %d)\n", #expr, __FILE__, __LINE__); \
      NAMD_die(msg); \
    } \
  } while (0)
#else
#define ASSERT(expr)
#endif 


// employ mixed precision
// (but allow easy change to all double precision for comparison)
typedef float Float;
typedef double Double;


namespace msm {

  ///////////////////////////////////////////////////////////////////////////
  //
  // Resizable Array class
  //
  ///////////////////////////////////////////////////////////////////////////

  template <class T> class Array;

  template <class T>
  void swap(Array<T>& s, Array<T>& t);

  template <class T>
  class Array {
    public:
      Array() : abuffer(0), alen(0), amax(0) { }
      Array(int n) : abuffer(0), alen(0), amax(0) { resize(n); }
      Array(const Array& a) : abuffer(0), alen(0), amax(0) { copy(a); }
      ~Array() { setmax(0); }
      Array& operator=(const Array& a) {
        if (this != &a) copy(a);  // don't allow self-assignment
        return(*this);
      }
      int len() const { return alen; }
      int max() const { return amax; }
      const T& operator[](int i) const {
#ifdef DEBUG_MSM
        return elem(i);
#else
        return abuffer[i];
#endif
      }
      const T& elem(int i) const {
        if (i < 0 || i >= alen) {
          char msg[100];
          snprintf(msg, sizeof(msg), "Array index:  alen=%d, i=%d\n", alen, i);
          NAMD_die(msg);
        }
        return abuffer[i];
      }
      T& operator[](int i) {
#ifdef DEBUG_MSM
        return elem(i);
#else
        return abuffer[i];
#endif
      }
      T& elem(int i) {
        if (i < 0 || i >= alen) {
          char msg[100];
          snprintf(msg, sizeof(msg), "Array index:  alen=%d, i=%d\n", alen, i);
          NAMD_die(msg);
        }
        return abuffer[i];
      }
      void append(const T& t) {
        if (alen==amax) setmax(2*amax+1);
        abuffer[alen++] = t;
      }
      void resize(int n) {
        if (n > amax) setmax(n);
        alen = n;
      }
      void setmax(int m);
      const T *buffer() const { return abuffer; }
      T *buffer() { return abuffer; }
      const T *buffer(int& n) const { n = alen; return abuffer; }
      T *buffer(int& n) { n = alen; return abuffer; }
      friend void swap<T>(Array&, Array&);
#ifdef DEBUG_MSM
      void print(const char *s=0) const {
        if (s) printf("PRINTING DATA FOR ARRAY \"%s\":\n", s);
        printf("abuffer=%p\n  alen=%d  amax=%d\n",
            (void *) abuffer, alen, amax);
      }
#endif
    protected:
      T *abuffer;
      int alen, amax;
      void copy(const Array& a);
  };

  template <class T>
  void Array<T>::setmax(int m) {
    if (m == amax) return;
    else if (m > 0) {
      T *newbuffer = new T[m];
      if ( ! newbuffer) {
        char msg[100];
        snprintf(msg, sizeof(msg), "Can't allocate %lu KB\n",
            (unsigned long)(m * sizeof(T) / 1024));
        NAMD_die(msg);
      }
      if (alen > m) alen = m;  // new buffer is shorter than old buffer
      for (int i = 0;  i < alen;  i++) {
        newbuffer[i] = abuffer[i];
      }
      delete[] abuffer;
      abuffer = newbuffer;
      amax = m;
    }
    else {  // consider m == 0
      delete[] abuffer;
      abuffer = 0;
      alen = 0;
      amax = 0;
    }
  }

  template <class T>
  void Array<T>::copy(const Array<T>& a) {
    setmax(a.amax);
    alen = a.alen;
    for (int i = 0;  i < alen;  i++) {
      abuffer[i] = a.abuffer[i];
    }
  }

  // swap arrays without duplicating memory buffer
  template <class T>
  void swap(Array<T>& s, Array<T>& t) {
    T *tmpbuffer = s.abuffer;  s.abuffer = t.abuffer;  t.abuffer = tmpbuffer;
    tmpbuffer = 0;
    int tmpn = s.alen;  s.alen = t.alen;  t.alen = tmpn;
    tmpn = s.amax;  s.amax = t.amax;  t.amax = tmpn;
    tmpn = s.astate;  s.astate = t.astate;  t.astate = tmpn;
  }


  ///////////////////////////////////////////////////////////////////////////
  //
  // Grid is 3D lattice of grid points with user-definable index ranges.
  //
  ///////////////////////////////////////////////////////////////////////////

  // 3-integer vector, used for indexing from a 3D grid
  struct Ivec {
    int i, j, k;
    Ivec(int n=0) : i(n), j(n), k(n) { }
    Ivec(int ni, int nj, int nk) : i(ni), j(nj), k(nk) { }
    int operator==(const Ivec& n) { return(i==n.i && j==n.j && k==n.k); }
#ifdef MSM_MIGRATION
    virtual void pup(PUP::er& p) {
      p|i, p|j, p|k;
    }
#endif
  };

  // index range for 3D lattice of grid points
  class IndexRange {
    public:
      IndexRange() : nlower(0), nextent(0) { }
      void set(int pia, int pni, int pja, int pnj, int pka, int pnk) {
        ASSERT(pni >= 0 && pnj >= 0 && pnk >= 0);
        nlower = Ivec(pia, pja, pka);
        nextent = Ivec(pni, pnj, pnk);
      }
      void setbounds(int pia, int pib, int pja, int pjb, int pka, int pkb) {
        set(pia, pib-pia+1, pja, pjb-pja+1, pka, pkb-pka+1);
      }
      int ia() const { return nlower.i; }
      int ib() const { return nlower.i + nextent.i - 1; }
      int ja() const { return nlower.j; }
      int jb() const { return nlower.j + nextent.j - 1; }
      int ka() const { return nlower.k; }
      int kb() const { return nlower.k + nextent.k - 1; }
      int ni() const { return nextent.i; }
      int nj() const { return nextent.j; }
      int nk() const { return nextent.k; }
      int nn() const { return nextent.i * nextent.j * nextent.k; }
      Ivec lower() const { return nlower; }
      Ivec extent() const { return nextent; }
      int operator<=(const IndexRange& n) {
        // true if this IndexRange fits inside n
        return ( ia() >= n.ia() && ib() <= n.ib() &&
                 ja() >= n.ja() && jb() <= n.jb() &&
                 ka() >= n.ka() && kb() <= n.kb() );
      }
#ifdef MSM_MIGRATION
      virtual void pup(PUP::er& p) {
        p|nlower, p|nextent;
      }
#endif
    protected:
      Ivec nlower;   // index for lowest corner of rectangular lattice
      Ivec nextent;  // extent of lattice along each dimension
  };

  // storage and indexing for 3D lattice of grid points
  // with fixed buffer storage no larger than size of block
  template <class T> class Grid;

  template <class T>
  class GridFixed : public IndexRange {
    friend class Grid<T>;
    public:
      GridFixed() { }
      void init(const IndexRange& n) {
        nlower = n.lower();
        nextent = n.extent();
        ASSERT(nextent.i * nextent.j * nextent.k <= MSM_MAX_BLOCK_VOLUME);
      }
      void set(int pia, int pni, int pja, int pnj, int pka, int pnk) {
        IndexRange::set(pia, pni, pja, pnj, pka, pnk);
        ASSERT(nextent.i * nextent.j * nextent.k <= MSM_MAX_BLOCK_VOLUME);
      }
      void setbounds(int pia, int pib, int pja, int pjb, int pka, int pkb) {
        IndexRange::setbounds(pia, pib, pja, pjb, pka, pkb);
        ASSERT(nextent.i * nextent.j * nextent.k <= MSM_MAX_BLOCK_VOLUME);
      }
      const T& operator()(int i, int j, int k) const {
#ifdef DEBUG_MSM
        return elem(i,j,k);
#else
        return gdata[flatindex(i,j,k)];
#endif
      }
      const T& operator()(const Ivec& n) const {
        return this->operator()(n.i, n.j, n.k);
      }
      const T& elem(int i, int j, int k) const {
        if (i<ia() || i>ib() || j<ja() || j>jb() || k<ka() || k>kb()) {
          char msg[200];
          snprintf(msg, sizeof(msg), "Grid indexing:\n"
              "ia=%d, ib=%d, i=%d\n"
              "ja=%d, jb=%d, j=%d\n"
              "ka=%d, kb=%d, k=%d\n",
              ia(), ib(), i, ja(), jb(), j, ka(), kb(), k);
          NAMD_die(msg);
        }
        return gdata[flatindex(i,j,k)];
      }
      T& operator()(int i, int j, int k) {
#ifdef DEBUG_MSM
        return elem(i,j,k);
#else
        return gdata[flatindex(i,j,k)];
#endif
      }
      T& operator()(const Ivec& n) {
        return this->operator()(n.i, n.j, n.k);
      }
      T& elem(int i, int j, int k) {
        if (i<ia() || i>ib() || j<ja() || j>jb() || k<ka() || k>kb()) {
          char msg[200];
          snprintf(msg, sizeof(msg), "Grid indexing:\n"
              "ia=%d, ib=%d, i=%d\n"
              "ja=%d, jb=%d, j=%d\n"
              "ka=%d, kb=%d, k=%d\n",
              ia(), ib(), i, ja(), jb(), j, ka(), kb(), k);
          NAMD_die(msg);
        }
        return gdata[flatindex(i,j,k)];
      }
      int flatindex(int i, int j, int k) const {
        return ((k-ka())*nj() + (j-ja()))*ni() + (i-ia());
      }
      const T *buffer() const { return gdata; }
      T *buffer() { return gdata; }

      // use to zero out grid
      void reset(const T& t) {
        int len = nn();
        for (int n = 0;  n < len;  n++) { gdata[n] = t; }
      }

      // use to modify the indexing by changing lower corner
      void updateLower(const Ivec& n) { nlower = n; }

      // accumulate another grid into this grid
      // the grid to be added must fit within this grid's index range
      GridFixed<T>& operator+=(const GridFixed<T>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        const T *gbuf = g.gdata.buffer();
        T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              buf[i + ijkoff] += gbuf[index];
            }
          }
        }
        return(*this);
      }

      // extract a subgrid from this grid
      // subgrid must fit within this grid's index range
      void extract(GridFixed<T>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        T *gbuf = g.gdata.buffer();
        const T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              gbuf[index] = buf[i + ijkoff];
            }
          }
        }
      }

    private:
      T gdata[MSM_MAX_BLOCK_VOLUME];
  };

  // storage and indexing for 3D lattice of grid points
  template <class T>
  class Grid : public IndexRange {
    public:
      Grid() { }
      void init(const IndexRange& n) {
        nlower = n.lower();
        nextent = n.extent();
        gdata.resize(nn());
      }
      void set(int pia, int pni, int pja, int pnj, int pka, int pnk) {
        IndexRange::set(pia, pni, pja, pnj, pka, pnk);
        gdata.resize(nn());
      }
      void setbounds(int pia, int pib, int pja, int pjb, int pka, int pkb) {
        IndexRange::setbounds(pia, pib, pja, pjb, pka, pkb);
        gdata.resize(nn());
      }
      void resize(int n) { // reserve space but don't set grid indexing
        gdata.resize(n);
      }
      const T& operator()(int i, int j, int k) const {
#ifdef DEBUG_MSM
        return elem(i,j,k);
#else
        return gdata[flatindex(i,j,k)];
#endif
      }
      const T& operator()(const Ivec& n) const {
        return this->operator()(n.i, n.j, n.k);
      }
      const T& elem(int i, int j, int k) const {
        if (i<ia() || i>ib() || j<ja() || j>jb() || k<ka() || k>kb()) {
          char msg[200];
          snprintf(msg, sizeof(msg), "Grid indexing:\n"
              "ia=%d, ib=%d, i=%d\n"
              "ja=%d, jb=%d, j=%d\n"
              "ka=%d, kb=%d, k=%d\n",
              ia(), ib(), i, ja(), jb(), j, ka(), kb(), k);
          NAMD_die(msg);
        }
        return gdata[flatindex(i,j,k)];
      }
      T& operator()(int i, int j, int k) {
#ifdef DEBUG_MSM
        return elem(i,j,k);
#else
        return gdata[flatindex(i,j,k)];
#endif
      }
      T& operator()(const Ivec& n) {
        return this->operator()(n.i, n.j, n.k);
      }
      T& elem(int i, int j, int k) {
        if (i<ia() || i>ib() || j<ja() || j>jb() || k<ka() || k>kb()) {
          char msg[200];
          snprintf(msg, sizeof(msg), "Grid indexing:\n"
              "ia=%d, ib=%d, i=%d\n"
              "ja=%d, jb=%d, j=%d\n"
              "ka=%d, kb=%d, k=%d\n",
              ia(), ib(), i, ja(), jb(), j, ka(), kb(), k);
          NAMD_die(msg);
        }
        return gdata[flatindex(i,j,k)];
      }
      int flatindex(int i, int j, int k) const {
        return ((k-ka())*nj() + (j-ja()))*ni() + (i-ia());
      }
      const Array<T>& data() const { return gdata; }
      Array<T>& data() { return gdata; }

      // use to zero out grid
      void reset(const T& t) {
        T *buf = gdata.buffer();
        int len = nn();
        for (int n = 0;  n < len;  n++) { buf[n] = t; }
      }

      // use to modify the indexing by changing lower corner
      void updateLower(const Ivec& n) { nlower = n; }

      // accumulate another grid into this grid
      // the grid to be added must fit within this grid's index range
      Grid<T>& operator+=(const Grid<T>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        const T *gbuf = g.gdata.buffer();
        T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              buf[i + ijkoff] += gbuf[index];
            }
          }
        }
        return(*this);
      }

      // extract a subgrid from this grid
      // subgrid must fit within this grid's index range
      void extract(Grid<T>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        T *gbuf = g.gdata.buffer();
        const T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              gbuf[index] = buf[i + ijkoff];
            }
          }
        }
      }

      // accumulate a fixed size grid into this grid
      // the grid to be added must fit within this grid's index range
      Grid<T>& operator+=(const GridFixed<T>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        const T *gbuf = g.buffer();
        T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              buf[i + ijkoff] += gbuf[index];
            }
          }
        }
        return(*this);
      }

      // extract a subgrid from this grid
      // subgrid must fit within this grid's index range
      void extract(GridFixed<T>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gni = g.nextent.i;
        int gnj = g.nextent.j;
        int gnk = g.nextent.k;
        int index = 0;
        int ni = nextent.i;
        int nij = nextent.i * nextent.j;
        int koff = (g.nlower.k - nlower.k) * nij
          + (g.nlower.j - nlower.j) * ni + (g.nlower.i - nlower.i);
        T *gbuf = g.buffer();
        const T *buf = gdata.buffer();
        for (int k = 0;  k < gnk;  k++) {
          int jkoff = k * nij + koff;
          for (int j = 0;  j < gnj;  j++) {
            int ijkoff = j * ni + jkoff;
            for (int i = 0;  i < gni;  i++, index++) {
              gbuf[index] = buf[i + ijkoff];
            }
          }
        }
      }

    private:
      Array<T> gdata;
  };


  ///////////////////////////////////////////////////////////////////////////
  //
  // Map object 
  //
  ///////////////////////////////////////////////////////////////////////////

  // index a block from the MSM grid hierarchy
  struct BlockIndex {
    int level;
    Ivec n;
    BlockIndex() : level(0), n(0) { }
    BlockIndex(int ll, const Ivec& nn) : level(ll), n(nn) { }
#ifdef MSM_MIGRATION
    virtual void pup(PUP::er& p) {
      p|level, p|n;
    }
#endif
  };

  // for uppermost levels of hierarchy
  // fold out image charges along periodic boundaries
  // to fill up desired block size
  struct FoldFactor {
    int active;   // is some numrep dimension > 1?
    Ivec numrep;  // number of replications along each dimension
    FoldFactor() : active(0), numrep(1) { }
    FoldFactor(int i, int j, int k) { set(i,j,k); }
    void set(int i, int j, int k) {
      if (i <= 0) i = 1;
      if (j <= 0) j = 1;
      if (k <= 0) k = 1;
      if (i > 1 || j > 1 || k > 1) active = 1;
      numrep = Ivec(i, j, k);
    }
  };

  // sending part of an extended grid calculation to another block
  struct BlockSend {
    BlockIndex nblock;       // relative block index
    IndexRange nrange;       // relative grid index range
    BlockIndex nblock_wrap;  // true block index
    IndexRange nrange_wrap;  // true grid index range
    void reset() {
      nblock = BlockIndex();
      nrange = IndexRange();
      nblock_wrap = BlockIndex();
      nrange_wrap = IndexRange();
    } // reset
#ifdef MSM_MIGRATION
    virtual void pup(PUP::er& p) {
      p|nblock, p|nrange, p|nblock_wrap, p|nrange_wrap;
    }
#endif
  };

  struct PatchSend {
    IndexRange nrange;         // true grid index range from my block
    IndexRange nrange_unwrap;  // relative grid index range for patch
    int patchID;               // patch ID
    void reset() {
      nrange = IndexRange();
      nrange_unwrap = IndexRange();
      patchID = -1;
    } // reset
  };

  // one PatchDiagram for each patch
  // maintain a Grid of PatchDiagram, indexed by patch ID
  struct PatchDiagram {
    IndexRange nrange;       // shows subset of MSM h-grid covering this patch
    Array<BlockSend> send;   // array of blocks to which this patch sends
    int numRecvs;            // number of blocks from which this patch receives
    void reset() {
      nrange = IndexRange();
      send.resize(0);
      numRecvs = 0;
    } // reset
  };

  // one BlockDiagram for each block of each level of each MSM grid
  // maintain a Grid of BlockDiagram for each level
  struct BlockDiagram {
    IndexRange nrange;            // subset of MSM grid for this block
    IndexRange nrangeCutoff;      // expanded subgrid for cutoff calculation
    IndexRange nrangeRestricted;  // (level+1) subgrid for restriction
    IndexRange nrangeProlongated; // (level-1) subgrid for prolongation
    Array<BlockSend> sendUp;      // send up charge to blocks on (level+1)
    Array<BlockSend> sendAcross;  // send across potential to blocks on (level)
    Array<int> indexGridCutoff;   // index of MsmGridCutoff chare to calculate
                                  // each charge -> potential block interaction
    Array<BlockSend> sendDown;    // send down potential to blocks on (level-1)
    Array<PatchSend> sendPatch;   // send my (level=0) potential block to patch
    int numRecvsCharge;           // number of expected receives of charge
    int numRecvsPotential;        // number of expected receives of potential

    void reset() {
      nrange = IndexRange();
      nrangeCutoff = IndexRange();
      nrangeRestricted = IndexRange();
      nrangeProlongated = IndexRange();
      sendUp.resize(0);
      sendAcross.resize(0);
      sendDown.resize(0);
      sendPatch.resize(0);
      numRecvsCharge = 0;
      numRecvsPotential = 0;
    } // reset
  };


  struct Map {
    Array<IndexRange> gridrange;  // dimensions for each MSM grid level

    Array<Grid<Float> > gc;     // grid constant weights for each level

    Array<PatchDiagram> patchList;
    Array<Grid<BlockDiagram> > blockLevel;

    int ispx, ispy, ispz;         // is periodic in x, y, z?

    Array<int> bsx, bsy, bsz;     // block size in x, y, z for each level

    Array<FoldFactor> foldfactor; // for uppermost grid levels
      // replicate periodic dimensions in order to fill up block size

    // clip index to grid level, using periodicity flags
    Ivec clipIndexToLevel(const Ivec& n, int level) const {
      ASSERT(level >= 0 && level < gridrange.len());
      Ivec pn(n);
      if ( ! ispx) {
        int a = gridrange[level].ia();
        int b = gridrange[level].ib();
        if (pn.i < a) pn.i = a;
        if (pn.i > b) pn.i = b;
      }
      if ( ! ispy) {
        int a = gridrange[level].ja();
        int b = gridrange[level].jb();
        if (pn.j < a) pn.j = a;
        if (pn.j > b) pn.j = b;
      }
      if ( ! ispz) {
        int a = gridrange[level].ka();
        int b = gridrange[level].kb();
        if (pn.k < a) pn.k = a;
        if (pn.k > b) pn.k = b;
      }
      return pn;
    }

    // determine relative (unwrapped) block index for the given grid index
    BlockIndex blockOfGridIndex(const Ivec& n, int level) const {
      ASSERT(level >= 0 && level < gridrange.len());
      BlockIndex bn;
      // we want floor((i - ia) / bsx), etc.
      // modify case i < ia to avoid integer division of negative numbers
      int d = n.i - gridrange[level].ia();
      bn.n.i = (d >= 0 ? d / bsx[level] : -((-d+bsx[level]-1) / bsx[level]));
      d = n.j - gridrange[level].ja();
      bn.n.j = (d >= 0 ? d / bsy[level] : -((-d+bsy[level]-1) / bsy[level]));
      d = n.k - gridrange[level].ka();
      bn.n.k = (d >= 0 ? d / bsz[level] : -((-d+bsz[level]-1) / bsz[level]));
      bn.level = level;
      return bn;
    }

    // determine relative (unwrapped) block index for the given grid index
    // for unfolded replication of image charges
    BlockIndex blockOfGridIndexFold(const Ivec& n, int level) const {
      ASSERT(level >= 0 && level < gridrange.len());
      BlockIndex bn;
      int bsi = foldfactor[level].numrep.i * bsx[level];
      int bsj = foldfactor[level].numrep.j * bsy[level];
      int bsk = foldfactor[level].numrep.k * bsz[level];
      // we want floor((i - ia) / bsx), etc.
      // modify case i < ia to avoid integer division of negative numbers
      int d = n.i - gridrange[level].ia();
      bn.n.i = (d >= 0 ? d / bsi : -((-d+bsi-1) / bsi));
      d = n.j - gridrange[level].ja();
      bn.n.j = (d >= 0 ? d / bsj : -((-d+bsj-1) / bsj));
      d = n.k - gridrange[level].ka();
      bn.n.k = (d >= 0 ? d / bsk : -((-d+bsk-1) / bsk));
      bn.level = level;
      return bn;
    }

    // find the natural index range of the given relative block number
    IndexRange indexRangeOfBlock(const BlockIndex& nb) const {
      ASSERT(nb.level >= 0 && nb.level < gridrange.len());
      IndexRange nr;
      int ia = nb.n.i * bsx[nb.level] + gridrange[nb.level].ia();
      int ja = nb.n.j * bsy[nb.level] + gridrange[nb.level].ja();
      int ka = nb.n.k * bsz[nb.level] + gridrange[nb.level].ka();
      nr.set(ia, bsx[nb.level], ja, bsy[nb.level], ka, bsz[nb.level]);
      return nr;
    }

    // find the natural index range of the given relative block number
    // for unfolded replication of image charges
    IndexRange indexRangeOfBlockFold(const BlockIndex& nb) const {
      ASSERT(nb.level >= 0 && nb.level < gridrange.len());
      int bsi = foldfactor[nb.level].numrep.i * bsx[nb.level];
      int bsj = foldfactor[nb.level].numrep.j * bsy[nb.level];
      int bsk = foldfactor[nb.level].numrep.k * bsz[nb.level];
      IndexRange nr;
      int ia = nb.n.i * bsi + gridrange[nb.level].ia();
      int ja = nb.n.j * bsj + gridrange[nb.level].ja();
      int ka = nb.n.k * bsk + gridrange[nb.level].ka();
      nr.set(ia, bsi, ja, bsj, ka, bsk);
      return nr;
    }

    // clip the natural block index range to not exceed the given index range
    IndexRange clipBlockToIndexRange(const BlockIndex& nb,
        const IndexRange& nrange) const {
      IndexRange nr = indexRangeOfBlock(nb);
      int nia = nrange.ia();
      int nib = nrange.ib();
      int nja = nrange.ja();
      int njb = nrange.jb();
      int nka = nrange.ka();
      int nkb = nrange.kb();
      int ia = nr.ia();
      if (ia < nia) ia = nia;
      int ib = nr.ib();
      if (ib > nib) ib = nib;
      int ja = nr.ja();
      if (ja < nja) ja = nja;
      int jb = nr.jb();
      if (jb > njb) jb = njb;
      int ka = nr.ka();
      if (ka < nka) ka = nka;
      int kb = nr.kb();
      if (kb > nkb) kb = nkb;
      nr.setbounds(ia, ib, ja, jb, ka, kb);
      return nr;
    }

    // clip the natural block index range to not exceed the given index range
    // for unfolded replication of image charges
    IndexRange clipBlockToIndexRangeFold(const BlockIndex& nb,
        const IndexRange& nrange) const {
      IndexRange nr = indexRangeOfBlockFold(nb);
      int nia = nrange.ia();
      int nib = nrange.ib();
      int nja = nrange.ja();
      int njb = nrange.jb();
      int nka = nrange.ka();
      int nkb = nrange.kb();
      int ia = nr.ia();
      if (ia < nia) ia = nia;
      int ib = nr.ib();
      if (ib > nib) ib = nib;
      int ja = nr.ja();
      if (ja < nja) ja = nja;
      int jb = nr.jb();
      if (jb > njb) jb = njb;
      int ka = nr.ka();
      if (ka < nka) ka = nka;
      int kb = nr.kb();
      if (kb > nkb) kb = nkb;
      nr.setbounds(ia, ib, ja, jb, ka, kb);
      return nr;
    }

    // set the nblock_wrap and nrange_wrap fields based on periodicity
    void wrapBlockSend(BlockSend& bs) const {
      BlockIndex nb = bs.nblock;
      IndexRange nr = bs.nrange;
      int level = bs.nblock.level;
      ASSERT(level >= 0 && level < blockLevel.len());
      int ni = blockLevel[level].ni();
      int nj = blockLevel[level].nj();
      int nk = blockLevel[level].nk();
      int di=0, dj=0, dk=0;
      if (ispx) {
        while (nb.n.i < 0) {
          nb.n.i += ni;
          di += ni * bsx[level];
        }
        while (nb.n.i >= ni) {
          nb.n.i -= ni;
          di -= ni * bsx[level];
        }
      }
      if (ispy) {
        while (nb.n.j < 0) {
          nb.n.j += nj;
          dj += nj * bsy[level];
        }
        while (nb.n.j >= nj) {
          nb.n.j -= nj;
          dj -= nj * bsy[level];
        }
      }
      if (ispz) {
        while (nb.n.k < 0) {
          nb.n.k += nk;
          dk += nk * bsz[level];
        }
        while (nb.n.k >= nk) {
          nb.n.k -= nk;
          dk -= nk * bsz[level];
        }
      }
      int ia = nr.ia();
      int ib = nr.ib();
      int ja = nr.ja();
      int jb = nr.jb();
      int ka = nr.ka();
      int kb = nr.kb();
      nr.setbounds(ia + di, ib + di, ja + dj, jb + dj, ka + dk, kb + dk);
      bs.nblock_wrap = nb;
      bs.nrange_wrap = nr;
    }

    // set the nblock_wrap and nrange_wrap fields based on periodicity
    // for unfolded replication of image charges
    void wrapBlockSendFold(BlockSend& bs) const {
      BlockIndex nb = bs.nblock;
      IndexRange nr = bs.nrange;
      int level = bs.nblock.level;
      ASSERT(level >= 0 && level < blockLevel.len());
      int foldi = foldfactor[level].numrep.i;
      int foldj = foldfactor[level].numrep.j;
      int foldk = foldfactor[level].numrep.k;
      int ni = blockLevel[level].ni();
      int nj = blockLevel[level].nj();
      int nk = blockLevel[level].nk();
      int bsi = foldi * bsx[level];
      int bsj = foldj * bsy[level];
      int bsk = foldk * bsz[level];
      int di=0, dj=0, dk=0;
      if (ispx) {
        while (nb.n.i < 0) {
          nb.n.i += ni;
          di += ni * bsi;
        }
        while (nb.n.i >= ni) {
          nb.n.i -= ni;
          di -= ni * bsi;
        }
      }
      if (ispy) {
        while (nb.n.j < 0) {
          nb.n.j += nj;
          dj += nj * bsj;
        }
        while (nb.n.j >= nj) {
          nb.n.j -= nj;
          dj -= nj * bsj;
        }
      }
      if (ispz) {
        while (nb.n.k < 0) {
          nb.n.k += nk;
          dk += nk * bsk;
        }
        while (nb.n.k >= nk) {
          nb.n.k -= nk;
          dk -= nk * bsk;
        }
      }
      int ia = nr.ia();
      int ib = nr.ib();
      int ja = nr.ja();
      int jb = nr.jb();
      int ka = nr.ka();
      int kb = nr.kb();
      nr.setbounds(ia + di, ib + di, ja + dj, jb + dj, ka + dk, kb + dk);
      bs.nblock_wrap = nb;
      bs.nrange_wrap = nr;
    }
  }; // Map


  struct AtomCoord {
    Position position;
    Real charge;
    int id;
  };

  typedef Array<AtomCoord> AtomCoordArray;
  typedef Array<Force> ForceArray;

  struct PatchData;
  typedef Array<PatchData> PatchDataArray;
  typedef Array<PatchData *> PatchPtrArray;

  struct BlockData;
  typedef Array<Grid<BlockData> > BlockDataGrids;

} // namespace msm

#endif // MSMMAP_H
