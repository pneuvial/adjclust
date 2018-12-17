//heap functions in C

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

//#define VERBOSE

// "row" indices in the chained list for each fusion (between a "left" cluster and a "right" cluster)
#define MINCL1 1  // min/leftmost position (in the original data) of the left cluster
#define MAXCL1 2  // max/rightmost position (in the original data) of the left cluster
#define MINCL2 3  // leftmost position (in the original data) of the right cluster
#define MAXCL2 4  // max/rightmost position (in the original data) of the right cluster
#define LAB1 5    // label of the left cluster (why call it "1" instead of "left"?)
#define LAB2 6    // label of the right cluster (why call it "2" instead of "right"?)
#define POSL 7    // position (in the chained list) of the left *neighbor* of the fusion (= a fusion!)
#define POSR 8    // position (in the chained list) of the right *neighbor* of the fusion (= a fusion!)
// NB: POSL and POSR are the elements that make it possible to navigate in the chained list
#define SII 9     // sum_{i,j in left cluster} sim(i,j) (why call it "II" instead of "left")
#define SJJ 10    // sum_{i,j in right cluster} sim(i,j)
#define SIJ 11    // sum_{i in left cluster, j in right cluster} sim(i,j)
#define VALID 12  // boolean indicator of the "validity" of a fusion (ie whether it involves cluster that still exist)

#define CHAIN(i,j) chainedL[(12*(j-1)+i-1)]
#define CUMR(i,j,p) rcCumRight[(p*(j-1)+i-1)]
#define CUML(i,j,p) rcCumLeft[(p*(j-1)+i-1)]
#define MERGE(i,j,k) merge[((*p-1)*(j-1)+i-1)]
#define TRW(i,j,k) traceW[((*p-1)*(j-1)+i-1)]

#define MIN(a,b) ((a) < (b) ? a : b)

//@export
SEXP percDown(SEXP Rpositions, SEXP Rdistances, SEXP Rl, SEXP Rpos){
    int mc, right, left;
    int *positions, *l, *pos;
    double *distances, tmp, val;
    Rpositions = PROTECT(coerceVector(Rpositions, INTSXP));
    positions = INTEGER(Rpositions);
    distances = REAL(Rdistances);
    l = INTEGER(Rl);
    pos = INTEGER(Rpos);
    *pos = *pos - 1;
    val = distances[positions[*pos]-1];
    
    while((2**pos+1) < *l) {
        //i=i+1;
        if ((2**pos+2) == *l){
            left = 2**pos+1;
            if (val > distances[positions[left]-1]) {  // positions in C start from 0 !!!
                // swap positions
                tmp = positions[*pos];
                positions[*pos] = positions[left];
                positions[left] = tmp;
                // update pos
                *pos = left;
            }
            else
                *pos = *l;
        }
        else {
            left = 2**pos+1;
            right = 2**pos+2;
            mc = right;
            if (distances[positions[left]-1] < distances[positions[right]-1]) // positions in C start from 0 !!!
                mc = left;
            
            if (val > distances[positions[mc]-1]) {  // positions in C start from 0 !!!
                // swap positions
                tmp = positions[*pos];
                positions[*pos] = positions[mc];
                positions[mc] = tmp;
                // update pos
                *pos = mc;
            }
            else
                *pos = *l;
        }
    }
    UNPROTECT (1) ;
    // return R_NilValue;
    return(Rpositions);
}

int* deleteMin_C(int *positions, double *distances, int l){
    int mc, right, left;
    int pos;
    double tmp, val;
    
    // pre-processing
    positions[0] = positions[l-1];
    l = l-1;
    pos = 1;
    
    pos = pos - 1;
    val = distances[positions[pos]-1];
    
    while((2*pos+1) < l) {
        // i=i+1;
        if ((2*pos+2) == l){
            left = 2*pos+1;
            if (val > distances[positions[left]-1]) {  // positions in C start from 0 !!!
                // swap positions
                tmp = positions[pos];
                positions[pos] = positions[left];
                positions[left] = tmp;
                // update pos
                pos = left;
            }
            else
                pos = l;
        }
        else {
            left = 2*pos+1;
            right = 2*pos+2;
            mc = right;
            if (distances[positions[left]-1] < distances[positions[right]-1]) // positions in C start from 0 !!!
                mc = left;
            
            if (val > distances[positions[mc]-1]) {  // positions in C start from 0 !!!
                // swap positions
                tmp = positions[pos];
                positions[pos] = positions[mc];
                positions[mc] = tmp;
                // update pos
                pos = mc;
            }
            else
                pos = l;
        }
    }
    // printf("%d\n", i);
    return(positions);
}

int* insertHeap_C(int *positions, double *distances, int l, int key){
    int pos;
    int parent;
    double tmp, val;
    
    // pre-processing
    positions[l] = key;
    pos = l+1;
    //*l = *l+1;
    
    parent = pos/2;
    val = distances[positions[pos-1]-1];
    
    while((pos-1) > 0){
        if ( distances[positions[parent-1]-1] > val){
            // i = i+1;
            // swap positions
            tmp = positions[parent-1];
            positions[parent-1] = positions[pos-1];
            positions[pos-1] = tmp;
            pos = pos/2;
            parent = pos/2;
        } else {
            pos=1;
        }
    }
    // printf("%d\n", i);
    
    return(positions);
}


// retrieve info for the right neighbor cluster of the fusion (via the right neighbor *fusion* of the fusion)
int* rightCluster_C(int posMin, double *chainedL){
    int *res = malloc(5*sizeof(int));  // vector of size 5
    if (CHAIN(POSR, posMin)<0) {
        res[0] = -1;
    } else  {
        res[0] = CHAIN(POSR, posMin);   
        res[1] = CHAIN(MINCL2, res[0]); 
        res[2] = CHAIN(MAXCL2, res[0]); 
        res[3] = CHAIN(LAB2, res[0]);   
        res[4] = CHAIN(POSR, res[0]);   
    } 
    return res;
}

// retrieve info for the left neighbor cluster of the fusion (via the left neighbor *fusion* of the fusion)
int* leftCluster_C(int posMin, double *chainedL){
    int *res = malloc(5*sizeof(int));  // vector of size 5
    if (CHAIN(POSL, posMin)<0){
        res[0] = -1;
    } else {  
        res[0] = CHAIN(POSL, posMin);   // index of the (left) *fusion* in the chained list of candidate fusions
        // the left neighbor cluster is the left element of the left fusion
        // info on the left neighbor cluster (of the left fusion):
        res[1] = CHAIN(MINCL1, res[0]); // min position of the 
        res[2] = CHAIN(MAXCL1, res[0]); // max position of the left neighbor cluster (of the left fusion)
        res[3] = CHAIN(LAB1, res[0]);   // label of the  left neighbor cluster (of the left fusion)
        // info on the left fusion (of the left fusion)
        res[4] = CHAIN(POSL, res[0]);   // position in the chained list of the left fusion (of the left fusion)
    }
    return res;
}

double rightPencil_C(int lim, int hLoc, int p, double *rcCumRight){
    double sumPen = CUMR(p, hLoc, p);
    if (lim<p){
        sumPen = sumPen - CUMR(p-lim, hLoc, p);
    }
    return(sumPen);
}

double leftPencil_C(int lim, int hLoc, int p, double *rcCumLeft){
    double sumPen = CUML(p, hLoc, p);
    if (lim>1){
        sumPen = sumPen - CUML(lim-1, hLoc, p);
    }
    return(sumPen);
}

double* distance_C(int mini, int maxi, int minj, int maxj, double *rcCumRight, double *rcCumLeft, int h, int p){
    double *res = malloc(sizeof(double)*4);
    double Sii, Sjj, Sij, D;
    int ni, nj, mIJ;
    
    ni = maxi - mini + 1;
    nj = maxj - minj + 1;
    mIJ = maxj - mini + 1;
    
    Sii = leftPencil_C(mini, MIN(h+1, ni), p, rcCumLeft) + 
        rightPencil_C(maxi, MIN(h+1, ni), p, rcCumRight) - 
        CUML(p, MIN(h+1, ni), p);
    
    Sjj = leftPencil_C(minj, MIN(h+1, nj), p, rcCumLeft) + 
        rightPencil_C(maxj, MIN(h+1, nj), p, rcCumRight) - 
        CUML(p, MIN(h+1, nj), p);
    
    Sij = (rightPencil_C(maxj, MIN(h+1, mIJ), p, rcCumRight) + 
        leftPencil_C(mini, MIN(h+1, mIJ), p, rcCumLeft) - 
        CUML(p, MIN(h+1, mIJ), p) - Sii - Sjj)/2;
    
    D = (float)ni*nj/(ni+nj)  * ( (float)1/(ni*ni)*Sii + (float)1/(nj*nj)*Sjj - (float)2/(ni*nj)*Sij ) ;
    
    res[0] = D;
    res[1] = Sii;
    res[2] = Sjj;
    res[3] = Sij;
    
#ifdef VERBOSE
    printf("mini: %d, ", mini);
    printf("maxi: %d, ", maxi);
    printf("minj: %d, ", minj);
    printf("maxj: %d\n", maxj);
    printf("Sii: %f, ", Sii);
    printf("Sjj: %f, ", Sjj);
    printf("Sij: %f\n", Sij);
    printf("D: %f\n", D);
#endif

    return res;
}

SEXP cWardHeaps(SEXP RrcCumRight, SEXP RrcCumLeft, SEXP Rh, SEXP Rp, SEXP RchainedL, SEXP Rpositions, SEXP Rdistances, SEXP RlHeap, SEXP Rmerge, SEXP Rgains, SEXP RtraceW){
    
    int *h, *p, *positions, *lHeap, *merge, *neiL, *neiR;
    int posMin, k;
    double *rcCumRight, *rcCumLeft, *distances, *chainedL, *gains, *traceW, *d1, *d2, *dLast, newDR, newDL, sumSdiag, snew, nii, njj, min_cl1, max_cl1, min_cl2, max_cl2;
    int jj, step, stepInv;
    
    Rpositions = PROTECT(coerceVector(Rpositions, INTSXP));
    Rmerge = PROTECT(coerceVector(Rmerge, INTSXP));  // matrix
    Rgains = PROTECT(coerceVector(Rgains, REALSXP));  // vector
    Rdistances = PROTECT(coerceVector(Rdistances, REALSXP));
    RchainedL = PROTECT(coerceVector(RchainedL, REALSXP));
    RtraceW = PROTECT(coerceVector(RtraceW, REALSXP));
    
    rcCumRight = REAL(RrcCumRight);
    rcCumLeft = REAL(RrcCumLeft);
    h = INTEGER(Rh);
    p = INTEGER(Rp);
    chainedL = REAL(RchainedL);
    positions = INTEGER(Rpositions);
    distances = REAL(Rdistances);
    lHeap = INTEGER(RlHeap);
    merge = INTEGER(Rmerge);
    gains = REAL(Rgains);
    traceW = REAL(RtraceW);
    
    k = *p - 1;
    
    jj = *p;
    sumSdiag = (float)0; // within cluster dispersion (WG)
    
    for ( step=1; step < (*p-1); step=step+1 ){
//    for ( step=1; step < 2; step=step+1 ){  // used as a hack to retrieve the heap status after the first merge...
        while(CHAIN(VALID, positions[0])==0){
            // printf("%f\n", distances[positions[0]-1]);
            // printf("%d\n", *lHeap);
            positions = deleteMin_C(positions, distances, *lHeap);
            *lHeap = *lHeap -1;
        }
        posMin = positions[0];   // index of the best fusion in the chained list
        // retrieve indices (in the set of elements to cluster) 
        // related to the left and right clusters of the best fusion:
        min_cl1 = CHAIN(MINCL1, posMin);  // leftmost element of the left cluster
        max_cl1 = CHAIN(MAXCL1, posMin);  // rightmost element of the left cluster
        min_cl2 = CHAIN(MINCL2, posMin);  // leftmost element of the right cluster
        max_cl2 = CHAIN(MAXCL2, posMin);  // rightmost element of the right cluster
        gains[step-1] = distances[posMin-1];  // current minimal value of the heap
        
        //remove head
        positions = deleteMin_C(positions, distances, *lHeap);
        *lHeap = *lHeap - 1;
        
        neiL = leftCluster_C(posMin, chainedL);  // 
        neiR = rightCluster_C(posMin, chainedL);
        
        if (min_cl1==1){ // left element of the fusion is the leftmost cluster
            d1 = distance_C(min_cl1, max_cl2, neiR[1], neiR[2], rcCumRight, rcCumLeft, *h, *p);
            newDR = d1[0];
            CHAIN(MINCL1, jj) = min_cl1;
            CHAIN(MAXCL1, jj) = max_cl2;
            CHAIN(MINCL2, jj) = neiR[1];
            CHAIN(MAXCL2, jj) = neiR[2];
            CHAIN(LAB1, jj) = step;
            CHAIN(LAB2, jj) = neiR[3];
            CHAIN(POSL, jj) = -1;
            CHAIN(POSR, jj) = neiR[4];
            CHAIN(SII, jj) = d1[1];
            CHAIN(SJJ, jj) = d1[2];
            CHAIN(SIJ, jj) = d1[3];
            CHAIN(VALID, jj) = 1;
            if (neiR[4]>0){
                CHAIN(POSL, neiR[4]) = jj;
            }
            CHAIN(VALID, posMin) = 0;
            CHAIN(VALID, neiR[0]) = 0;
            distances[jj-1] = newDR;
            positions = insertHeap_C(positions, distances, *lHeap, jj);
            *lHeap = *lHeap + 1;
            jj = jj+1;
        }
        else if (max_cl2==*p){ // right element of the fusion is the rightmost cluster
            d2 = distance_C(neiL[1], neiL[2], min_cl1, max_cl2, rcCumRight, rcCumLeft, *h, *p);
            newDL = d2[0] ;
            CHAIN(MINCL1, jj) = neiL[1];
            CHAIN(MAXCL1, jj) = neiL[2];
            CHAIN(MINCL2, jj) = min_cl1;
            CHAIN(MAXCL2, jj) = max_cl2;
            CHAIN(LAB1, jj) = neiL[3];
            CHAIN(LAB2, jj) = step;
            CHAIN(POSL, jj) = neiL[4];
            CHAIN(POSR, jj) = -1;
            CHAIN(SII, jj) = d2[1];
            CHAIN(SJJ, jj) = d2[2];
            CHAIN(SIJ, jj) = d2[3];
            CHAIN(VALID, jj) = 1;
            if (neiL[4]>0){
                CHAIN(POSR, neiL[4]) = jj;
            }
            CHAIN(VALID, posMin) = 0;
            CHAIN(VALID, neiL[0]) = 0;
            distances[jj-1] = newDL;
            positions = insertHeap_C(positions, distances, *lHeap, jj);
            *lHeap = *lHeap + 1;
            jj = jj+1;
        }
        else {  // fusion does not involve leftmost or rightmost cluster
            d2 = distance_C(min_cl1, max_cl2, neiR[1], neiR[2], rcCumRight, rcCumLeft, *h, *p);
            d1 = distance_C(neiL[1], neiL[2], min_cl1, max_cl2, rcCumRight, rcCumLeft, *h, *p);
            newDR = d2[0];
            newDL = d1[0];
            CHAIN(MINCL1, jj) = neiL[1];
            CHAIN(MAXCL1, jj) = neiL[2];
            CHAIN(MINCL2, jj) = min_cl1;
            CHAIN(MAXCL2, jj) = max_cl2;
            CHAIN(LAB1, jj) = neiL[3];  // name of left cluster of the new fusion (between left of best fusion and best fusion) = name of the left fusion of the best fusion
            CHAIN(LAB2, jj) = step;     // name of right cluster of the new fusion (between left of best fusion and best fusion) = name of the cluster created by the best fusion, ie a new name
            CHAIN(POSL, jj) = neiL[4];  // index of fusion in chained list
            CHAIN(POSR, jj) = jj+1;     // index of fusion in chained list
            CHAIN(SII, jj) = d1[1];
            CHAIN(SJJ, jj) = d1[2];
            CHAIN(SIJ, jj) = d1[3];
            CHAIN(VALID, jj) = 1;
            CHAIN(MINCL1, jj+1) = min_cl1;
            CHAIN(MAXCL1, jj+1) = max_cl2;
            CHAIN(MINCL2, jj+1) = neiR[1];
            CHAIN(MAXCL2, jj+1) = neiR[2];
            CHAIN(LAB1, jj+1) = step;
            CHAIN(LAB2, jj+1) = neiR[3];
            CHAIN(POSL, jj+1) = jj;
            CHAIN(POSR, jj+1) = neiR[4];
            CHAIN(SII, jj+1) = d2[1];
            CHAIN(SJJ, jj+1) = d2[2];
            CHAIN(SIJ, jj+1) = d2[3];
            CHAIN(VALID, jj+1) = 1;
            if (neiL[4]>0){
                CHAIN(POSR, neiL[4]) = jj;
            }
            if (neiR[4]>0){
                CHAIN(POSL, neiR[4]) = jj+1;
            }
            CHAIN(VALID, posMin) = 0;
            CHAIN(VALID, neiL[0]) = 0;
            CHAIN(VALID, neiR[0]) = 0;
            distances[jj-1] = newDL;
            positions = insertHeap_C(positions, distances, *lHeap, jj);
            *lHeap = *lHeap + 1;
            distances[jj] = newDR;
            positions = insertHeap_C(positions, distances, *lHeap, jj+1);
            *lHeap = *lHeap + 1;
            jj = jj+2;
        }
        // update of sumSdiag
        nii = max_cl1 - min_cl1 + 1;
        njj = max_cl2 - min_cl2 + 1;
        snew = CHAIN(SII, posMin) + CHAIN(SJJ, posMin) + (float)2*CHAIN(SIJ, posMin);
        sumSdiag = CHAIN(SII, posMin)/nii + CHAIN(SJJ, posMin)/njj - snew/(nii+njj);
        
        stepInv = *p - step;
        MERGE(step, 1, k) = CHAIN(LAB1, posMin);
        MERGE(step, 2, k) = CHAIN(LAB2, posMin);
        TRW(stepInv, 1, k) = stepInv;
        TRW(stepInv, 2, k) = (float)sumSdiag;
    }
    
    // merging the last two classes
    step = *p-1;
    stepInv = 1;
    posMin = jj-1;
    min_cl1 = CHAIN(MINCL1, posMin);
    max_cl1 = CHAIN(MAXCL1, posMin);
    min_cl2 = CHAIN(MINCL2, posMin);
    max_cl2 = CHAIN(MAXCL2, posMin);
    dLast = distance_C(min_cl1, max_cl1, min_cl2, max_cl2, rcCumRight, rcCumLeft, *h, *p);
    // update of sumSdiag
    nii = max_cl1 - min_cl1 + 1;
    njj = max_cl2 - min_cl2 + 1;
    snew = CHAIN(SII, posMin) + CHAIN(SJJ, posMin) + (float)2*CHAIN(SIJ, posMin);
    sumSdiag = CHAIN(SII, posMin)/nii + CHAIN(SJJ, posMin)/njj - snew/(nii+njj);
    
    gains[step-1] = dLast[0];
    MERGE(step, 1, k) = CHAIN(LAB1, jj-1);
    MERGE(step, 2, k) = CHAIN(LAB2, jj-1);
    TRW(stepInv, 1, k) = stepInv;
    TRW(stepInv, 2, k) = (float)sumSdiag;
    
    UNPROTECT(6);
    return(Rmerge);
}
