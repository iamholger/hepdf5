#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include <highfive/H5Easy.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

using HighFive::Chunking;
using HighFive::AtomicType;
using HighFive::Deflate;
using HighFive::File;
using HighFive::DataSpace;
using HighFive::DataSetCreateProps;

#include "HepMC3/ReaderAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"
#include "HepMC3/GenParticle.h"

using namespace HepMC3;
using namespace std;

namespace HepDF5
{
    vector<string> split(const string& s, char delimiter)
    {
        vector<string> tokens;
        string token;
        istringstream tokenStream(s);
        while (getline(tokenStream, token, delimiter))
        {
            tokens.push_back(token);
        }
        return tokens;
    }

    vector<string> weightnames_from_event(GenEvent const & evt)
    {
        return evt.weight_names();
    }

    vector<string> weightnames_from_ascii_file(string const & fname)
    {
        ReaderAsciiHepMC2 reader(fname);
        if (reader.failed()) 
        {
        }

        GenEvent evt;
        reader.read_event(evt);
        vector<string> wnames = weightnames_from_event(evt);
        reader.close();
        return wnames; 
    }

    int hepmc_version_from_file(string const & fname)
    {

        int rc(0);
        ifstream infile(fname);
        string sLine;
        if (infile.good())
        {
            while (getline(infile, sLine))
            {
                if (sLine.rfind("HepMC::Version", 0) == 0)
                {
                    auto tokens = split(sLine, ' ');
                    auto nums = split(tokens[1], '.');
                    rc=atoi(nums[0].c_str());
                    break;
                }
            }
        }
        infile.close();
        return rc;
    }

    void create_datasets(HighFive::File  & file, vector<string> wnames = {}, int compression=4, size_t chunksize=32)
    {
        DataSetCreateProps props;
        props.add(Deflate(compression));
        props.add(Chunking(std::vector<hsize_t>{chunksize}));
        
        file.createGroup("particle");
        file.createGroup("vertex");
        file.createGroup("event");
        
        vector<string> dsname_dbl = 
        { 
            "particle/px", "particle/py", "particle/pz", "particle/E",
            "particle/m", "vertex/x", "vertex/y", "vertex/z", "vertex/t",
        };
        for (auto dsname : dsname_dbl)
        {
            file.createDataSet(dsname, DataSpace( {0}, {DataSpace::UNLIMITED}), AtomicType<double>(), props);
        }

        vector<string> dsname_int =
        { 
            "particle/status", "particle/pid", "particle/vtx_prod",
            "particle/vtx_end", "vertex/id", "vertex/status"
        };
        for (auto dsname : dsname_int)
        {
            file.createDataSet(dsname, DataSpace( {0}, {DataSpace::UNLIMITED}), AtomicType<int>(), props);
        }

        vector<string> dsname_szt =
        { 
            "event/start_p", "event/start_v", "event/npart", "event/nvert"
        };
        for (auto dsname : dsname_szt)
        {
            file.createDataSet(dsname, DataSpace( {0}, {DataSpace::UNLIMITED}), AtomicType<size_t>(), props);
        }

        H5Easy::dump(file, "weight_names", wnames);
        DataSetCreateProps wprops;
        wprops.add(Deflate(compression));
        wprops.add(Chunking(std::vector<hsize_t>{32,wnames.size()}));
        file.createDataSet("event/weight",
            DataSpace( 
              {0, wnames.size()}, 
              {DataSpace::UNLIMITED,
              wnames.size()}), AtomicType<double>(), wprops);
    }

    void write_particles(HighFive::File  & file, vector<GenParticlePtr> const & parts)
    {
      size_t const nparts = parts.size();
      size_t current = file.getDataSet("particle/status").getSpace().getDimensions()[0];
      size_t after = current + nparts;

      
      for (auto k : file.getGroup("particle").listObjectNames())
      {
        file.getDataSet("particle/" + k).resize({after});
      }

      vector<int> status, pid, vtx_end, vtx_prod;
      vector<double> e, m, px, py, pz;
      status.reserve(nparts);
      pid.reserve(nparts);
      vtx_end.reserve(nparts);
      vtx_prod.reserve(nparts);
      for (auto p : parts)
      {
        status.push_back(p->status());
        pid.push_back(p->pid());
        if (p->production_vertex() != NULL)
        {
          vtx_prod.push_back(p->production_vertex()->id());
        }
        else
        {
          vtx_prod.push_back(-1);
        }
        
        if (p->end_vertex() != NULL)
        {
          vtx_end.push_back(p->end_vertex()->id());
        }
        else
        {
          vtx_end.push_back(-1);
        }

        e.push_back(p->momentum().e());
        m.push_back(p->momentum().m());
        px.push_back(p->momentum().px());
        py.push_back(p->momentum().py());
        pz.push_back(p->momentum().pz());
      }
      file.getDataSet("particle/status")  .select({current}, {nparts}).write(status);
      file.getDataSet("particle/pid")     .select({current}, {nparts}).write(pid);
      file.getDataSet("particle/vtx_end") .select({current}, {nparts}).write(vtx_end);
      file.getDataSet("particle/vtx_prod").select({current}, {nparts}).write(vtx_prod);
      file.getDataSet("particle/E").select({current}, {nparts}).write(e);
      file.getDataSet("particle/m").select({current}, {nparts}).write(m);
      file.getDataSet("particle/px").select({current}, {nparts}).write(px);
      file.getDataSet("particle/py").select({current}, {nparts}).write(py);
      file.getDataSet("particle/pz").select({current}, {nparts}).write(pz);

    } 

    void write_vertices(HighFive::File  & file, vector<GenVertexPtr> const & verts)
    {
      size_t const nvtx = verts.size();
      size_t current = file.getDataSet("vertex/status").getSpace().getDimensions()[0];
      size_t after = current + nvtx;
      
      for (auto k : file.getGroup("vertex").listObjectNames())
      {
        file.getDataSet("vertex/" + k).resize({after});
      }

      vector<int> status, vid;
      status.reserve(nvtx);
      vid.reserve(nvtx);

      vector<double> vx, vy, vz, vt;
      vx.reserve(nvtx);
      vy.reserve(nvtx);
      vz.reserve(nvtx);
      vt.reserve(nvtx);

      for (auto v : verts)
      {
        status.push_back(v->status());
        vid.push_back(v->id());
        vx.push_back(v->position().x());
        vy.push_back(v->position().y());
        vz.push_back(v->position().z());
        vt.push_back(v->position().t());
      }
      file.getDataSet("vertex/status")  .select({current}, {nvtx}).write(status);
      file.getDataSet("vertex/id")     .select({current}, {nvtx}).write(vid);
      file.getDataSet("vertex/x")     .select({current}, {nvtx}).write(vx);
      file.getDataSet("vertex/y")     .select({current}, {nvtx}).write(vy);
      file.getDataSet("vertex/z")     .select({current}, {nvtx}).write(vz);
      file.getDataSet("vertex/t")     .select({current}, {nvtx}).write(vt);
    }

    size_t load_part(File const & file, string const & path, size_t idx) 
    {
       auto dataset = file.getDataSet(path);
       vector<size_t> data;
       dataset.select({idx}, {1}).read(data);
       return data[0];
    }

    void write_event(HighFive::File  & file, vector< vector<double> > WGT, vector<size_t> NP, vector<size_t> NV)
    {
      size_t const nevt = NV.size();
      size_t current = file.getDataSet("event/npart").getSpace().getDimensions()[0];
      size_t after = current + nevt;
      
      for (auto k : file.getGroup("event").listObjectNames())
      {
        if (k=="weight") continue;
        file.getDataSet("event/" + k).resize({after});
      }
      file.getDataSet("event/weight").resize({after, WGT[0].size()});

      size_t pstart(0), vstart(0);

      if (current>0)
      {
        pstart = load_part(file, "event/start_p", current-1) + load_part(file, "event/npart", current-1);
        vstart = load_part(file, "event/start_v", current-1) + load_part(file, "event/nvert", current-1);
      }

      file.getDataSet("event/npart")  .select({current}, {nevt}).write(NP);
      file.getDataSet("event/nvert")  .select({current}, {nevt}).write(NV);

      vector<size_t> cs_part(nevt);
      cs_part[0] = 0;
      partial_sum(NP.begin(), NP.end()-1, cs_part.begin()+1, plus<size_t>());
      transform(cs_part.begin(), cs_part.end(), cs_part.begin(), bind2nd( plus<size_t>(), pstart));

      vector<size_t> cs_vtx(nevt);
      cs_vtx[0] = 0;
      partial_sum(NV.begin(), NV.end()-1, cs_vtx.begin()+1, plus<size_t>());
      transform(cs_vtx.begin(), cs_vtx.end(), cs_vtx.begin(), bind2nd( plus<size_t>(), vstart));

      file.getDataSet("event/start_p")  .select({current}, {nevt}).write(cs_part);
      file.getDataSet("event/start_v")  .select({current}, {nevt}).write(cs_vtx);
     
      file.getDataSet("event/weight")  .select({current, 0}, {nevt, WGT[0].size()}).write(WGT);
    } 

    void test(string const & fname, HighFive::File  & fout, size_t write_every=100000, size_t msg_every=100)
    {

        int iev(0);
        size_t npart(0), nvtx(0);

        vector<size_t> nV, nP;

        // TODO reserve resize etc
        vector<GenParticlePtr> wP;
        vector<GenVertexPtr>   wV;
        vector< vector<double> >   wW;

        ReaderAsciiHepMC2 reader(fname);
        while (!reader.failed()) 
        {
            GenEvent evt;
            reader.read_event(evt);
            if (reader.failed())  {
                printf("End of file reached. Exit.\n");
                break;
            }

            auto P = evt.particles();
            npart=P.size();
            for (auto p : P) wP.push_back(p);
            nP.push_back(npart);

            auto V = evt.vertices();
            nvtx = V.size();
            for (auto v : V) wV.push_back(v);
            nV.push_back(nvtx);
            
            auto W = evt.weights();
            wW.push_back(W);

            evt.clear();

            if (wP.size() > write_every)
            {
              write_particles(fout, wP);
              write_vertices(fout, wV);
              write_event(fout, wW, nP, nV);
              wP.clear();
              wV.clear();
              wW.clear();
              nV.clear();
              nP.clear();
            }
            
            iev++;
            if ((iev+1)%msg_every==0)
            {
              cerr << "Events processed: " << iev+1 << "\r"; 
            }

        }
        if (wP.size() > 0)
        {
            write_particles(fout, wP);
            write_vertices(fout, wV);
            write_event(fout, wW, nP, nV);
        }

    }

}

int main(int argc, const char* argv[]) 
{
    auto version = HepDF5::hepmc_version_from_file(argv[1]);
    if (version==0)
    {
        cerr << "Unable to determine hepmc version from file " << argv[1] << " exiting.";
        exit(1);
    }
    else if (version==2)
    {
        cerr << "HepMC version 2 detected on input file\n";
        vector<string> const WN = HepDF5::weightnames_from_ascii_file(argv[1]);
        H5Easy::File outfile("myexample.h5", H5Easy::File::ReadWrite | H5Easy::File::Create| H5Easy::File::Truncate);
        HepDF5::create_datasets(outfile, WN, atoi(argv[3]), atoi(argv[4]));
        HepDF5::test(argv[1], outfile, atoi(argv[2]));
    }
    else if (version ==3)
    {
        cerr << "HepMC version 3 detected on input file\n";
    }
    else
    {
        cerr << "Unsupported HepMC version " << version << " detected on input file, exiting\n";
        exit(1);
    }

    cerr << "I am done.\n";
    return 0;
}

