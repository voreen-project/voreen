
namespace voreen {
   struct VesselGraph;
}

struct NetmetsResult {
   float fnr;
   float fpr;
};

extern "C" {
   NetmetsResult netmets_compare_networks(const voreen::VesselGraph& g1, const voreen::VesselGraph& g2);
}
