#pragma once

namespace hepdf5
{

  std::vector<std::tuple<size_t, size_t> > distribute(size_t nevents, size_t nread, int nranks, int rank)
  {
    std::vector<std::tuple<size_t, size_t> > strips;

    size_t nfirst=0;
    size_t ntot=0;
    size_t readsize = nread;
    if (nranks==1)
    {
      while (ntot < nevents)
      {
        if (nevents-ntot < readsize) readsize=nevents-ntot;
        {
        }
        strips.push_back({ntot+nfirst, readsize});
        ntot+=readsize;
      }

    }
    else
    {
      size_t rankread = nevents/nranks;
      nfirst = rank*rankread;
      while (ntot < rankread)
      {
        if (rankread-ntot < readsize) readsize=rankread-ntot;
        {
        }
        strips.push_back({ntot+nfirst, readsize});
        ntot+=readsize;
      }
    }

    if (rank+1 == nranks)
    {
      std::get<1>(strips.back()) = nevents - std::get<0>(strips.back());
    }

    return strips;
  }

  struct Keeper
  {
      Keeper(size_t nevents, size_t nread, int nranks, int rank) 
          :
              strips(distribute(nevents, nread, nranks, rank)), 
              n_strip(strips.size())
      {
          auto [os, rs] = strips[0];
          n_ev = rs;
          i_ev=0;
          i_strip=0;
      }

      std::unique_ptr<int> next() 
      {
          if (i_ev == n_ev)
          {
              if (i_strip==n_strip-1) return NULL;
              auto [offset, rs] = strips[++i_strip];
              // fillBuffes(offset,rs)
              n_ev+=rs;
          }
          return std::make_unique<int>(i_ev++);
      }

      std::vector<std::tuple<size_t, size_t> > const  strips;
      size_t const n_strip;
      size_t n_ev;
      size_t i_ev;
      size_t i_strip;
  };
}
