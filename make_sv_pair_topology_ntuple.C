#include <TFile.h>
#include <TTree.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

bool isHeavyFlavorStatus(int status) {
  return status > 1 && status % 100 == 0;
}

struct SvInfo {
  int svIndex = -1;
  int assignedCode = 0;
  int svtxNtrk = 0;
  float svtxdl = -999.0;
  float svtxdls = -999.0;
  float svtxdl2d = -999.0;
  float svtxdls2d = -999.0;
  float svtxm = -999.0;
  float svtxmcorr = -999.0;
  float svtxpt = -999.0;
  float svtxnormchi2 = -999.0;
  int nTracksUsed = 0;
  int nHeavyFlavorTracks = 0;
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double sumPt = 0.0;
  double maxAbsIp3dSig = 0.0;
  std::map<int, double> hfPtByStatus;
};

struct OutRow {
  Long64_t entry = -1;
  int jetIndex = -1;
  float jtpt = 0.0;
  float jteta = 0.0;
  float uPartB = 0.0;
  int truthCode1 = 0;
  int truthCode2 = 0;
  int svCode1 = 0;
  int svCode2 = 0;
  int topology = 0;  // 1 = 1-1, 2 = 2-0
  int is11 = 0;
  int is20 = 0;

  int svIndex1 = -1;
  int svIndex2 = -1;
  float svpt1 = 0.0;
  float svpt2 = 0.0;
  float svm1 = 0.0;
  float svm2 = 0.0;
  float svmcorr1 = 0.0;
  float svmcorr2 = 0.0;
  float svdls1 = 0.0;
  float svdls2 = 0.0;
  float svdls2d1 = 0.0;
  float svdls2d2 = 0.0;
  float svdl1 = 0.0;
  float svdl2 = 0.0;
  float svdl2d1 = 0.0;
  float svdl2d2 = 0.0;
  float svnormchi2_1 = 0.0;
  float svnormchi2_2 = 0.0;
  int svNtrk1 = 0;
  int svNtrk2 = 0;
  int svTrackCount1 = 0;
  int svTrackCount2 = 0;
  float maxAbsIp3dSig1 = 0.0;
  float maxAbsIp3dSig2 = 0.0;

  float svEta1 = 0.0;
  float svEta2 = 0.0;
  float svPhi1 = 0.0;
  float svPhi2 = 0.0;
  float deltaEta = 0.0;
  float deltaPhi = 0.0;
  float deltaR = 0.0;
  float ptRatio = 0.0;
  float massSum = 0.0;
  float massRatio = 0.0;
  float dlsDiff = 0.0;
  float dlsAbsDiff = 0.0;
  float dls2dDiff = 0.0;
  float dls2dAbsDiff = 0.0;
};

int assignSvStatusByPt(const SvInfo& sv) {
  int bestStatus = 0;
  double bestPt = -1.0;
  for (const auto& item : sv.hfPtByStatus) {
    if (item.second > bestPt) {
      bestPt = item.second;
      bestStatus = item.first;
    }
  }
  return bestStatus;
}

double etaFromP(double px, double py, double pz) {
  const double p = std::sqrt(px * px + py * py + pz * pz);
  if (p <= std::abs(pz)) return pz >= 0 ? 999.0 : -999.0;
  return 0.5 * std::log((p + pz) / (p - pz));
}

double phiFromP(double px, double py) {
  return std::atan2(py, px);
}

void addTrackToSv(SvInfo& sv, float pt, float eta, float phi, int status,
                  float ip3dSig = 0.0) {
  const double px = pt * std::cos(phi);
  const double py = pt * std::sin(phi);
  const double pz = pt * std::sinh(eta);
  sv.px += px;
  sv.py += py;
  sv.pz += pz;
  sv.sumPt += pt;
  ++sv.nTracksUsed;
  sv.maxAbsIp3dSig = std::max(sv.maxAbsIp3dSig, static_cast<double>(std::abs(ip3dSig)));
  if (isHeavyFlavorStatus(status)) {
    sv.hfPtByStatus[status] += pt;
    ++sv.nHeavyFlavorTracks;
  }
}

void fillSvFields(const SvInfo& sv, bool first, OutRow& row) {
  const float eta = etaFromP(sv.px, sv.py, sv.pz);
  const float phi = phiFromP(sv.px, sv.py);
  if (first) {
    row.svIndex1 = sv.svIndex;
    row.svCode1 = sv.assignedCode;
    row.svpt1 = sv.svtxpt;
    row.svm1 = sv.svtxm;
    row.svmcorr1 = sv.svtxmcorr;
    row.svdls1 = sv.svtxdls;
    row.svdls2d1 = sv.svtxdls2d;
    row.svdl1 = sv.svtxdl;
    row.svdl2d1 = sv.svtxdl2d;
    row.svnormchi2_1 = sv.svtxnormchi2;
    row.svNtrk1 = sv.svtxNtrk;
    row.svTrackCount1 = sv.nTracksUsed;
    row.maxAbsIp3dSig1 = sv.maxAbsIp3dSig;
    row.svEta1 = eta;
    row.svPhi1 = phi;
  } else {
    row.svIndex2 = sv.svIndex;
    row.svCode2 = sv.assignedCode;
    row.svpt2 = sv.svtxpt;
    row.svm2 = sv.svtxm;
    row.svmcorr2 = sv.svtxmcorr;
    row.svdls2 = sv.svtxdls;
    row.svdls2d2 = sv.svtxdls2d;
    row.svdl2 = sv.svtxdl;
    row.svdl2d2 = sv.svtxdl2d;
    row.svnormchi2_2 = sv.svtxnormchi2;
    row.svNtrk2 = sv.svtxNtrk;
    row.svTrackCount2 = sv.nTracksUsed;
    row.maxAbsIp3dSig2 = sv.maxAbsIp3dSig;
    row.svEta2 = eta;
    row.svPhi2 = phi;
  }
}

void makeBranches(TTree& tree, OutRow& row) {
  tree.Branch("entry", &row.entry, "entry/L");
  tree.Branch("jetIndex", &row.jetIndex, "jetIndex/I");
  tree.Branch("jtpt", &row.jtpt, "jtpt/F");
  tree.Branch("jteta", &row.jteta, "jteta/F");
  tree.Branch("uPartB", &row.uPartB, "uPartB/F");
  tree.Branch("truthCode1", &row.truthCode1, "truthCode1/I");
  tree.Branch("truthCode2", &row.truthCode2, "truthCode2/I");
  tree.Branch("svCode1", &row.svCode1, "svCode1/I");
  tree.Branch("svCode2", &row.svCode2, "svCode2/I");
  tree.Branch("topology", &row.topology, "topology/I");
  tree.Branch("is11", &row.is11, "is11/I");
  tree.Branch("is20", &row.is20, "is20/I");

  tree.Branch("svIndex1", &row.svIndex1, "svIndex1/I");
  tree.Branch("svIndex2", &row.svIndex2, "svIndex2/I");
  tree.Branch("svpt1", &row.svpt1, "svpt1/F");
  tree.Branch("svpt2", &row.svpt2, "svpt2/F");
  tree.Branch("svm1", &row.svm1, "svm1/F");
  tree.Branch("svm2", &row.svm2, "svm2/F");
  tree.Branch("svmcorr1", &row.svmcorr1, "svmcorr1/F");
  tree.Branch("svmcorr2", &row.svmcorr2, "svmcorr2/F");
  tree.Branch("svdls1", &row.svdls1, "svdls1/F");
  tree.Branch("svdls2", &row.svdls2, "svdls2/F");
  tree.Branch("svdls2d1", &row.svdls2d1, "svdls2d1/F");
  tree.Branch("svdls2d2", &row.svdls2d2, "svdls2d2/F");
  tree.Branch("svdl1", &row.svdl1, "svdl1/F");
  tree.Branch("svdl2", &row.svdl2, "svdl2/F");
  tree.Branch("svdl2d1", &row.svdl2d1, "svdl2d1/F");
  tree.Branch("svdl2d2", &row.svdl2d2, "svdl2d2/F");
  tree.Branch("svnormchi2_1", &row.svnormchi2_1, "svnormchi2_1/F");
  tree.Branch("svnormchi2_2", &row.svnormchi2_2, "svnormchi2_2/F");
  tree.Branch("svNtrk1", &row.svNtrk1, "svNtrk1/I");
  tree.Branch("svNtrk2", &row.svNtrk2, "svNtrk2/I");
  tree.Branch("svTrackCount1", &row.svTrackCount1, "svTrackCount1/I");
  tree.Branch("svTrackCount2", &row.svTrackCount2, "svTrackCount2/I");
  tree.Branch("maxAbsIp3dSig1", &row.maxAbsIp3dSig1, "maxAbsIp3dSig1/F");
  tree.Branch("maxAbsIp3dSig2", &row.maxAbsIp3dSig2, "maxAbsIp3dSig2/F");

  tree.Branch("svEta1", &row.svEta1, "svEta1/F");
  tree.Branch("svEta2", &row.svEta2, "svEta2/F");
  tree.Branch("svPhi1", &row.svPhi1, "svPhi1/F");
  tree.Branch("svPhi2", &row.svPhi2, "svPhi2/F");
  tree.Branch("deltaEta", &row.deltaEta, "deltaEta/F");
  tree.Branch("deltaPhi", &row.deltaPhi, "deltaPhi/F");
  tree.Branch("deltaR", &row.deltaR, "deltaR/F");
  tree.Branch("ptRatio", &row.ptRatio, "ptRatio/F");
  tree.Branch("massSum", &row.massSum, "massSum/F");
  tree.Branch("massRatio", &row.massRatio, "massRatio/F");
  tree.Branch("dlsDiff", &row.dlsDiff, "dlsDiff/F");
  tree.Branch("dlsAbsDiff", &row.dlsAbsDiff, "dlsAbsDiff/F");
  tree.Branch("dls2dDiff", &row.dls2dDiff, "dls2dDiff/F");
  tree.Branch("dls2dAbsDiff", &row.dls2dAbsDiff, "dls2dAbsDiff/F");
}

}  // namespace

void make_sv_pair_topology_ntuple(
    const char* input =
        "/Users/mnguyen/Downloads/merged_block_0000_Pythia8_recalJP_20260617.root",
    const char* treeName = "ak4PFJetAnalyzer/t",
    const char* output = "sv_pair_topology_pt80_140.root",
    Long64_t maxEntries = -1,
    double jetPtMin = 80.0,
    double jetPtMax = 140.0,
    double jetAbsEtaMax = 2.0,
    double minUParTB = -1.0) {
  TFile* file = TFile::Open(input);
  if (!file || file->IsZombie()) {
    std::cerr << "Could not open input file: " << input << "\n";
    return;
  }

  TTree* tree = nullptr;
  file->GetObject(treeName, tree);
  if (!tree) {
    std::cerr << "Could not find tree: " << treeName << "\n";
    return;
  }

  int nref = 0;
  float jtpt[1000];
  float jteta[1000];
  int jtNbHad[1000];
  float discrUParTProbB[1000];
  float discrUParTProbBB[1000];
  float discrUParTProbLeptonicB[1000];

  int nfullB = 0;
  int fullBJetId[1000];
  int fullBSta[1000];

  int nsvtx = 0;
  int svtxJetId[10000];
  int svtxNtrk[10000];
  float svtxdl[10000];
  float svtxdls[10000];
  float svtxdl2d[10000];
  float svtxdls2d[10000];
  float svtxm[10000];
  float svtxmcorr[10000];
  float svtxpt[10000];
  float svtxnormchi2[10000];

  int ntrk = 0;
  int trkJetId[20000];
  int trkSvtxId[20000];
  float trkPt[20000];
  float trkEta[20000];
  float trkPhi[20000];
  float trkIp3dSig[20000];
  int trkMatchSta[20000];

  int ntrkInSvtxNotInJet = 0;
  int trkInSvtxNotInJetSvId[20000];
  float trkInSvtxNotInJetPt[20000];
  float trkInSvtxNotInJetEta[20000];
  float trkInSvtxNotInJetPhi[20000];
  int trkInSvtxNotInJetMatchSta[20000];

  tree->SetBranchStatus("*", 0);
  const char* branches[] = {
      "nref", "jtpt", "jteta", "jtNbHad",
      "discr_unifiedParticleTransformer_probb",
      "discr_unifiedParticleTransformer_probbb",
      "discr_unifiedParticleTransformer_problepb",
      "nfullB", "fullBJetId", "fullBSta",
      "nsvtx", "svtxJetId", "svtxNtrk", "svtxdl", "svtxdls", "svtxdl2d",
      "svtxdls2d", "svtxm", "svtxmcorr", "svtxpt", "svtxnormchi2",
      "ntrk", "trkJetId", "trkSvtxId", "trkPt", "trkEta", "trkPhi",
      "trkIp3dSig", "trkMatchSta",
      "ntrkInSvtxNotInJet", "trkInSvtxNotInJetSvId", "trkInSvtxNotInJetPt",
      "trkInSvtxNotInJetEta", "trkInSvtxNotInJetPhi",
      "trkInSvtxNotInJetMatchSta"};
  for (const char* branch : branches) tree->SetBranchStatus(branch, 1);
  tree->SetCacheSize(256 * 1024 * 1024);
  tree->AddBranchToCache("*", true);

  tree->SetBranchAddress("nref", &nref);
  tree->SetBranchAddress("jtpt", jtpt);
  tree->SetBranchAddress("jteta", jteta);
  tree->SetBranchAddress("jtNbHad", jtNbHad);
  tree->SetBranchAddress("discr_unifiedParticleTransformer_probb", discrUParTProbB);
  tree->SetBranchAddress("discr_unifiedParticleTransformer_probbb", discrUParTProbBB);
  tree->SetBranchAddress("discr_unifiedParticleTransformer_problepb",
                         discrUParTProbLeptonicB);
  tree->SetBranchAddress("nfullB", &nfullB);
  tree->SetBranchAddress("fullBJetId", fullBJetId);
  tree->SetBranchAddress("fullBSta", fullBSta);
  tree->SetBranchAddress("nsvtx", &nsvtx);
  tree->SetBranchAddress("svtxJetId", svtxJetId);
  tree->SetBranchAddress("svtxNtrk", svtxNtrk);
  tree->SetBranchAddress("svtxdl", svtxdl);
  tree->SetBranchAddress("svtxdls", svtxdls);
  tree->SetBranchAddress("svtxdl2d", svtxdl2d);
  tree->SetBranchAddress("svtxdls2d", svtxdls2d);
  tree->SetBranchAddress("svtxm", svtxm);
  tree->SetBranchAddress("svtxmcorr", svtxmcorr);
  tree->SetBranchAddress("svtxpt", svtxpt);
  tree->SetBranchAddress("svtxnormchi2", svtxnormchi2);
  tree->SetBranchAddress("ntrk", &ntrk);
  tree->SetBranchAddress("trkJetId", trkJetId);
  tree->SetBranchAddress("trkSvtxId", trkSvtxId);
  tree->SetBranchAddress("trkPt", trkPt);
  tree->SetBranchAddress("trkEta", trkEta);
  tree->SetBranchAddress("trkPhi", trkPhi);
  tree->SetBranchAddress("trkIp3dSig", trkIp3dSig);
  tree->SetBranchAddress("trkMatchSta", trkMatchSta);
  tree->SetBranchAddress("ntrkInSvtxNotInJet", &ntrkInSvtxNotInJet);
  tree->SetBranchAddress("trkInSvtxNotInJetSvId", trkInSvtxNotInJetSvId);
  tree->SetBranchAddress("trkInSvtxNotInJetPt", trkInSvtxNotInJetPt);
  tree->SetBranchAddress("trkInSvtxNotInJetEta", trkInSvtxNotInJetEta);
  tree->SetBranchAddress("trkInSvtxNotInJetPhi", trkInSvtxNotInJetPhi);
  tree->SetBranchAddress("trkInSvtxNotInJetMatchSta", trkInSvtxNotInJetMatchSta);

  TFile* outFile = TFile::Open(output, "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "Could not create output file: " << output << "\n";
    return;
  }
  TTree outTree("svPairTopology", "Clean two-SV topology diagnostic");
  OutRow row;
  makeBranches(outTree, row);

  long long selectedJets = 0;
  long long exactlyTwoTruth = 0;
  long long exactlyTwoRecoSv = 0;
  long long bothAssignedTruth = 0;
  long long n11 = 0;
  long long n20 = 0;
  long long skippedNonTruthAssignment = 0;

  const Long64_t nEntries =
      maxEntries >= 0 ? std::min(maxEntries, tree->GetEntries()) : tree->GetEntries();

  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    if (entry > 0 && entry % 1000000 == 0) {
      std::cerr << "Processed " << entry << " / " << nEntries << " entries\n";
    }
    tree->GetEntry(entry);

    for (int ijet = 0; ijet < nref; ++ijet) {
      if (!(jtpt[ijet] > jetPtMin && jtpt[ijet] < jetPtMax)) continue;
      if (!(std::abs(jteta[ijet]) < jetAbsEtaMax)) continue;
      if (jtNbHad[ijet] != 2) continue;
      const double uPartB =
          discrUParTProbB[ijet] + discrUParTProbBB[ijet] + discrUParTProbLeptonicB[ijet];
      if (!(uPartB > minUParTB)) continue;
      ++selectedJets;

      std::set<int> truthCodes;
      for (int ib = 0; ib < nfullB; ++ib) {
        if (fullBJetId[ib] != ijet) continue;
        if (isHeavyFlavorStatus(fullBSta[ib])) truthCodes.insert(fullBSta[ib]);
      }
      if (truthCodes.size() != 2) continue;
      ++exactlyTwoTruth;

      std::map<int, SvInfo> svInfos;
      for (int isv = 0; isv < nsvtx; ++isv) {
        if (svtxJetId[isv] != ijet) continue;
        SvInfo& sv = svInfos[isv];
        sv.svIndex = isv;
        sv.svtxNtrk = svtxNtrk[isv];
        sv.svtxdl = svtxdl[isv];
        sv.svtxdls = svtxdls[isv];
        sv.svtxdl2d = svtxdl2d[isv];
        sv.svtxdls2d = svtxdls2d[isv];
        sv.svtxm = svtxm[isv];
        sv.svtxmcorr = svtxmcorr[isv];
        sv.svtxpt = svtxpt[isv];
        sv.svtxnormchi2 = svtxnormchi2[isv];
      }
      if (svInfos.size() != 2) continue;
      ++exactlyTwoRecoSv;

      for (int itrk = 0; itrk < ntrk; ++itrk) {
        if (trkJetId[itrk] != ijet) continue;
        if (!svInfos.count(trkSvtxId[itrk])) continue;
        addTrackToSv(svInfos[trkSvtxId[itrk]], trkPt[itrk], trkEta[itrk],
                     trkPhi[itrk], trkMatchSta[itrk], trkIp3dSig[itrk]);
      }

      for (int itrk = 0; itrk < ntrkInSvtxNotInJet; ++itrk) {
        const int svId = trkInSvtxNotInJetSvId[itrk];
        if (!svInfos.count(svId)) continue;
        addTrackToSv(svInfos[svId], trkInSvtxNotInJetPt[itrk],
                     trkInSvtxNotInJetEta[itrk], trkInSvtxNotInJetPhi[itrk],
                     trkInSvtxNotInJetMatchSta[itrk], 0.0);
      }

      std::vector<SvInfo> svs;
      for (auto& item : svInfos) {
        item.second.assignedCode = assignSvStatusByPt(item.second);
        svs.push_back(item.second);
      }
      if (svs[0].assignedCode == 0 || svs[1].assignedCode == 0) continue;
      if (!truthCodes.count(svs[0].assignedCode) || !truthCodes.count(svs[1].assignedCode)) {
        ++skippedNonTruthAssignment;
        continue;
      }
      ++bothAssignedTruth;

      std::sort(svs.begin(), svs.end(),
                [](const SvInfo& lhs, const SvInfo& rhs) {
                  if (lhs.svtxpt != rhs.svtxpt) return lhs.svtxpt > rhs.svtxpt;
                  return lhs.svIndex < rhs.svIndex;
                });

      row = OutRow();
      row.entry = entry;
      row.jetIndex = ijet;
      row.jtpt = jtpt[ijet];
      row.jteta = jteta[ijet];
      row.uPartB = uPartB;
      auto truthIt = truthCodes.begin();
      row.truthCode1 = *truthIt++;
      row.truthCode2 = *truthIt;

      fillSvFields(svs[0], true, row);
      fillSvFields(svs[1], false, row);

      row.is11 = svs[0].assignedCode != svs[1].assignedCode;
      row.is20 = svs[0].assignedCode == svs[1].assignedCode;
      row.topology = row.is11 ? 1 : 2;
      if (row.is11) ++n11;
      if (row.is20) ++n20;

      row.deltaEta = row.svEta1 - row.svEta2;
      row.deltaPhi = TVector2::Phi_mpi_pi(row.svPhi1 - row.svPhi2);
      row.deltaR = std::sqrt(row.deltaEta * row.deltaEta + row.deltaPhi * row.deltaPhi);
      const float minPt = std::min(row.svpt1, row.svpt2);
      const float maxPt = std::max(row.svpt1, row.svpt2);
      row.ptRatio = maxPt > 0.0 ? minPt / maxPt : 0.0;
      row.massSum = row.svm1 + row.svm2;
      const float minMass = std::min(row.svm1, row.svm2);
      const float maxMass = std::max(row.svm1, row.svm2);
      row.massRatio = maxMass > 0.0 ? minMass / maxMass : 0.0;
      row.dlsDiff = row.svdls1 - row.svdls2;
      row.dlsAbsDiff = std::abs(row.dlsDiff);
      row.dls2dDiff = row.svdls2d1 - row.svdls2d2;
      row.dls2dAbsDiff = std::abs(row.dls2dDiff);

      outTree.Fill();
    }
  }

  outFile->cd();
  outTree.Write();
  outFile->Close();

  std::cout << "Input: " << input << "\n";
  std::cout << "Output: " << output << "\n";
  std::cout << "Processed entries: " << nEntries << "\n";
  std::cout << "Jet cuts: " << jetPtMin << " < jtpt < " << jetPtMax
            << ", abs(jteta) < " << jetAbsEtaMax << ", jtNbHad == 2"
            << ", UParT b score > " << minUParTB << "\n";
  std::cout << "Selected jtNbHad==2 jets: " << selectedJets << "\n";
  std::cout << "Exactly two fullB truth codes: " << exactlyTwoTruth << "\n";
  std::cout << "Exactly two reco SVs: " << exactlyTwoRecoSv << "\n";
  std::cout << "Both SVs assigned to truth HF codes: " << bothAssignedTruth << "\n";
  std::cout << "Skipped non-truth HF assignments: " << skippedNonTruthAssignment << "\n";
  std::cout << "1-1 rows: " << n11 << "\n";
  std::cout << "2-0 rows: " << n20 << "\n";
}
