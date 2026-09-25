// Bin-by-bin comparison of two ROOT files, used to prove that re-running a job
// reproduces it exactly. ROOT files are never byte-identical (UUIDs, timestamps),
// so compare the contents instead.
//   root -l -b -q 'compare_root_contents.C("a.root","b.root")'
// Prints IDENTICAL / DIFFERENT and the first few differing bins.
void compare_root_contents(TString fa, TString fb)
{
    TFile *A = TFile::Open(fa), *B = TFile::Open(fb);
    if (!A || A->IsZombie() || !B || B->IsZombie()) { printf("cannot open inputs\n"); return; }

    int nhist = 0, ndiff = 0, nmissing = 0;
    double maxdiff = 0.;
    TIter next(A->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *) next())) {
        TObject *oa = key->ReadObj();
        TH1 *ha = dynamic_cast<TH1 *>(oa);
        if (!ha) continue;
        nhist++;
        TH1 *hb = dynamic_cast<TH1 *>(B->Get(key->GetName()));
        if (!hb) { printf("  MISSING in b: %s\n", key->GetName()); nmissing++; continue; }
        const int nb = ha->GetNcells();
        if (nb != hb->GetNcells()) { printf("  BINNING differs: %s\n", key->GetName()); ndiff++; continue; }
        int bad = 0;
        for (int i = 0; i < nb; ++i) {
            const double d = std::fabs(ha->GetBinContent(i) - hb->GetBinContent(i));
            const double e = std::fabs(ha->GetBinError(i)   - hb->GetBinError(i));
            if (d > 0. || e > 0.) { bad++; maxdiff = std::max(maxdiff, std::max(d, e)); }
        }
        if (bad) {
            ndiff++;
            if (ndiff <= 5)
                printf("  DIFFERS: %-45s %d/%d bins, entries %g vs %g\n",
                       key->GetName(), bad, nb, ha->GetEntries(), hb->GetEntries());
        }
    }
    printf("%s\n  vs %s\n", fa.Data(), fb.Data());
    printf("  %d histograms, %d differ, %d missing, max |delta| = %g  ->  %s\n",
           nhist, ndiff, nmissing, maxdiff,
           (ndiff == 0 && nmissing == 0) ? "IDENTICAL" : "DIFFERENT");
    A->Close(); B->Close();
}
