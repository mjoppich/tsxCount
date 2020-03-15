//
// Created by mjopp on 11.06.2016.
//

#ifndef TSXCOUNT_BIJECTIVEKMAPPING_H
#define TSXCOUNT_BIJECTIVEKMAPPING_H

#include "IBijectiveFunction.h"
#include <cstdlib>
#include <ctime>
#include <cstring>

class BijectiveKMapping;

class BijectiveMappingException: public std::exception
{
public:
    BijectiveMappingException(BijectiveKMapping* pMapping)
            : std::exception()
    {
        m_pMapping = pMapping;
    }


    virtual const char* what() const throw()
    {
        return this->getExceptText();
    }

    virtual const char* getExceptText() const = 0;

protected:

    BijectiveKMapping* m_pMapping = NULL;

};



class BijectiveMappingNoInverseException : BijectiveMappingException
{
public:

    BijectiveMappingNoInverseException(BijectiveKMapping* pMapping)
    : BijectiveMappingException(pMapping)
    {

    }

    virtual const char* getExceptText() const;
};

class BijectiveMapping1DeterminantException : BijectiveMappingException
{
public:

    BijectiveMapping1DeterminantException(BijectiveKMapping* pMapping)
            : BijectiveMappingException(pMapping)
    {

    }

    virtual const char* getExceptText() const;
};
/**
 * \class BijectiveKMapping
 *
 *
 * \brief An implementation of the
 *
 * This class implements a bijective function \see{IBijectiveFunction}
 * This class is used as a mapping from a kmer -> key value
 *
 * \author Markus Joppich
 *
 */
class BijectiveKMapping : public IBijectiveFunction {

public:

    BijectiveKMapping(uint32_t iK, MemoryPool<FIELDTYPE>* pPool)
    :m_iK(iK), m_matN(2*iK), m_pPool(pPool)
    {
        srand(time(NULL));

        m_pA = createMatrix(m_matN);
        m_pAi = createMatrix(m_matN);

        if (!inverse(m_pA, m_pAi, m_matN))
        {
            throw BijectiveMappingNoInverseException(this);
        }

        //std::cout << this->printMatrix(m_pA) << std::endl;

        m_pArows = this->matrixToRows(m_pA);
        m_pAirows = this->matrixToRows(m_pAi);

    }

    int8_t* getOneVector()
    {

        int8_t* pReturn = (int8_t*) calloc( sizeof(int8_t), m_matN );

        for (uint32_t i = 0; i < m_matN; ++i)
            pReturn[i] = 1;

        return pReturn;

    }

    bool equalVec(int8_t* pVec1, int8_t* pVec2)
    {
        for (uint32_t i = 0; i < m_matN; ++i)
        {
            if (pVec1[i] != pVec2[i])
                return false;
        }

        return true;
    }

    int mod (int a, int b)
    {
        if(b < 0) //you can check for b == 0 separately and do what you want
            return mod(a, -b);
        int ret = a % b;
        if(ret < 0)
            ret+=b;

        return ret;
    }

    /**
     * \brief applies a bijective function (in binary) to iKmer
     * \value iKmer the kmer bits to apply this function to
     * \return the return key
     */
    virtual UBigInt apply(UBigInt& iKmer){

        return applyto(iKmer, m_pArows);
    }
    /**
     * \brief applies a bijective function (in binary) to a key
     * \value iKey the key bits to apply this function to
     * \return the originating kmer
     */
    virtual UBigInt inv_apply(UBigInt& iKey)
    {
        return applyto(iKey, m_pAirows);
    }



    uint32_t getK() {
        return m_iK;
    }

    std::string printMatrix()
    {
        return this->printMatrix(m_pA);
    }

    std::string printInvMatrix()
    {
        return this->printMatrix(m_pAi);
    }

    template<class T>
    std::string printMatrix(T* A)
    {
        std::stringstream oSS;

        for (uint32_t i=0; i<m_matN; i++)
        {
            for (uint32_t j=0; j<m_matN; j++)
                oSS << std::to_string(A[pos(i,j)]) << " ";
            oSS << std::endl;
        }

        return oSS.str();
    }

    template<class T>
    std::string printVector(T* A)
    {
        std::stringstream oSS;

        oSS << "[ ";

        for (int i=0; i<m_matN; i++)
            oSS << std::to_string(A[i]) << " ";

        oSS << "]" << std::endl;

        return oSS.str();
    }

protected:

    UBigInt applyto(UBigInt& oValue, UBigInt* pRows)
    {
        UBigInt oReturn(m_matN, true, this->m_pPool);

        for (uint32_t i = 0; i < m_matN; ++i)
        {

            UBigInt oTmp(pRows[i], this->m_pPool);
            oTmp.bitAnd(oValue);
            //oRes = pRows[i] & oValue; // was ^ for whatever reason? (XOR)

            //std::cerr << pRows[i].to_string() << std::endl;
            //std::cerr << oValue.to_string() << std::endl;
            //std::cerr << oRes.to_string() << std::endl;

            uint32_t iBitSum = oTmp.sumBits();

            oReturn.setBit(m_matN-1-i, iBitSum % 2);


        }

        return oReturn;
    }

    UBigInt* matrixToRows( int8_t* pMatrix) {

        UBigInt* pReturn = (UBigInt*) ::malloc(sizeof(UBigInt) * 2* m_iK);//new UBigInt[m_iK * 2];

        for (uint32_t i = 0; i < m_iK * 2; ++i)
        {
            UBigInt newElem(m_iK*2, true, m_pPool);
            newElem.transferSelf(pReturn+i);

            //pReturn[i].resize(m_iK*2);
            //memcpy(pReturn + i, &oRet, sizeof(UBigInt));

            //pReturn[i] = oRet;

            for (uint32_t j = 0; j < m_iK * 2; ++j)
            {

                uint8_t iValue = pMatrix[pos(i,j)];

                pReturn[i].setBit(2*m_iK - 1 -j, iValue);

            }

            //std::cout << pReturn[i].to_string() << std::endl;


        }

        return pReturn;
    }

    void applyRingMat(int8_t* pA)
    {
        this->applyRing(pA, m_matN, m_matN);
    }

    void applyRingVec(int8_t* pA)
    {
        this->applyRing(pA, 1, m_matN);
    }

    void applyRing(int8_t* pA, int lenI, int lenJ)
    {
        for (int i = 0; i < lenI; ++i)
        {
            for (int j = 0; j < lenJ; ++j)
            {
                pA[pos(i,j)] = mod(pA[pos(i,j)], 2);
            }
        }
    }

    int8_t* getZeroMatrix()
    {
        return (int8_t*) calloc(sizeof(int8_t), m_matN*m_matN);
    }

    void getRandomMatrix(int8_t* pMat, size_t n)
    {

        srand (time(NULL));

        for (uint32_t i = 0; i < m_matN; ++i)
        {
            pMat[pos(i,i)] = 1;

            for (uint32_t j = i+1; j < m_matN; ++j)
            {

                int8_t iValue = rand() % 2;
                pMat[pos(i,j)] = iValue;

            }

        }

    }

    template <class T>
    int8_t* getVectorMul(T* pMatA, int8_t* pVec)
    {

        int8_t* pReturn = (int8_t*) calloc(sizeof(int8_t), m_matN);

        for (uint32_t i = 0; i < m_matN; ++i)
        {

            int8_t iValue = 0;

            for (uint32_t k = 0; k < m_matN; ++k)
            {

                iValue += pMatA[pos(i,k)] * pVec[k];

            }

            pReturn[i] = iValue;

        }

        return pReturn;

    }

    int8_t* getMatrixMul(int8_t* pMatA, int8_t* pMatB, size_t n)
    {

        int8_t* pReturn = getZeroMatrix();

        for (uint32_t i = 0; i < n; ++i)
        {

            for (uint32_t j = 0; j < n; ++j)
            {

                int8_t iSum = 0;

                for (uint32_t k = 0; k < n; ++k)
                {

                    iSum += pMatA[pos(i,k)] * pMatB[pos(k,j)];

                }

                pReturn[pos(i,j)] = iSum;

            }

        }

        return pReturn;
    }

    int8_t* createVec(size_t n)
    {
        return (int8_t*) calloc( sizeof(int8_t), n);
    }

    int8_t iabs(int8_t val)
    {

        if (val >= 0)
            return val;

        return -val;
    }

    bool decomposeMatrixRecipe(int8_t* a, int8_t* pivot, size_t n)
    {

        int8_t* vv = createVec(n);
        int8_t d = 1;

        // should not do scaling ;)

        int8_t big = 0;

        for (size_t i = 0; i < n; ++i)
        {
            pivot[i] = i;
        }

        for (size_t i = 0; i < n; ++i)
        {
            big = 0;

            for (size_t j = 1; j < n; ++j)
            {
                if (a[npos(i,j,n)] > big)
                    big = a[npos(i,j,n)];

            }

            if (big == 0)
            {
                std::cerr << "error in factors: singularity" << std::endl;

                free(vv);
                return false;
            }

            vv[i] = 1 / big;
        }

        for (size_t i = 0; i < n; ++i)
        {

            if (vv[i] != 1)
            {
                std::cerr << "error in factors: factor detected: " << std::to_string(vv[i]) << std::endl;

                free(vv);
                return false;
            }
        }

        //printVector(vv, n);

        int8_t imax = 0;

        for (size_t j = 0; j < n; ++j)
        {

            // from 0 to j
            for (size_t i = 0; i<j; ++i)
            {

                int8_t sum = a[npos(i,j,n)];
                for (size_t k = 0; k < i; ++k)
                {
                    sum -= a[npos(i,k,n)] * a[npos(k,j,n)];
                }
                a[npos(i,j,n)] = sum;

            }


            // from j to n
            big = 0;
            imax = j;
            for (size_t i = j; i < n; ++i)
            {
                int8_t sum = a[npos(i,j,n)];

                for (size_t k = 0; k < j; ++k)
                {
                    sum -= a[npos(i,k,n)] * a[npos(k,j,n)];
                }
                a[npos(i,j,n)] = sum;


                if (vv[i] * iabs(sum) > big)
                {
                    big = vv[i];
                    imax = i;
                }
            }


            // row pivotizing
            if (j != imax)
            {
                //std::cout << std::to_string(j) << std::endl;

                for (size_t k = 0; k < n; ++k)
                {
                    int8_t dum = a[npos(imax, k, n)];
                    a[npos(imax, k, n)] = a[npos(j,k, n)];
                    a[npos(j, k, n)] = dum;
                }

                d = -d;
                vv[imax] = vv[j];

                size_t tmp = pivot[j];
                pivot[j]=pivot[imax];
                pivot[imax] = tmp;

            }


            if (a[npos(j,j,n)] == 0)
            {
                //printVector(pivot, n);

                //std::cout << "jj == 0 " << std::to_string(j) << std::endl;
                free(vv);
                return false;
            }

            if (j != n)
            {

                int8_t dum = 1 / a[npos(j,j,n)];

                for (size_t i = j+1; i < n; ++i)
                {
                    a[npos(i,j,n)] *= dum;
                }

            }

        }

        free(vv);

        return true;

    }

    void lubksb(int8_t* mat, size_t n, int8_t* pivot, int8_t* vec)
    {

        size_t ip;
        //ii = 0;
        ip = 0;

        int8_t sum = 0;

        for (size_t i = 0; i < n; ++i)
        {

            ip = pivot[i];
            sum = vec[ip];
            vec[ip] = vec[i];

            //vec[i] = sum;
            if (i > 0)
            {
                for (size_t j = 0; j < i; ++j)
                {
                    //std::cout << "( " << std::to_string(i) << " , " << std::to_string(j) << " )" << std::endl;
                    sum -= mat[npos(i,j,n)] * vec[j];
                }
            }


            //std::cerr << "in i: " << std::to_string(sum) << std::endl;

            vec[i]=sum;

        }

        //std::cout << "2.3.6 ";
        //printVector(vec, n);

        for (int64_t i = n-1; i >= 0; --i)
        {

            sum = vec[i];

            for (size_t j = i+1; j < n; ++j)
                sum -= mat[npos(i,j,n)]*vec[j];

            if (mat[npos(i,i,n)] == 0)
                std::cerr << "div null " << std::to_string(i) << std::endl;

            vec[i] = sum / mat[npos(i,i,n)];

        }

        //std::cout << "2.3.7 ";
        //printVector(vec, n);

    }


    int8_t* invertComb(int8_t* comb, int8_t* pivot, size_t n)
    {

        int8_t* inv = createMatrix(n);
        int8_t* col = createVec(n);

        for (size_t j = 0; j < n; ++j)
        {

            for(size_t i = 0; i < n; ++i)
            {
                col[i] = 0;
            }

            col[j] = 1;

            lubksb(comb, n, pivot, col);

            for (size_t i = 0; i < n; ++i)
            {
                inv[npos(i,j,n)] = col[i];
            }

        }

        return inv;

    }

    int8_t triangularDeterminant(int8_t* pmat, size_t n)
    {

        int8_t det=1;

        for (size_t i = 0; i < n; ++i)
            det *= pmat[npos(i,i,n)];

        return det;

    }

    int8_t* matrixMul(int8_t* pL, int8_t* pR, size_t n)
    {
        int8_t* res = (int8_t*) calloc( sizeof(int8_t), n*n);

        for (size_t i = 0; i < n; ++i)
        {

            for (size_t j = 0; j < n; ++j)
            {

                int8_t sum = 0;

                for (size_t k = 0; k < n; ++k)
                {

                    sum += pL[npos(i,k,n)] * pR[npos(k,j,n)];

                }

                res[npos(i,j,n)] = sum;

            }

        }

        return res;
    }

    bool inverse(int8_t* pA, int8_t* pAi, size_t n)
    {
        int8_t* pMat;
        int8_t* LU;
        int8_t* pivot;

        int8_t iDeterminant = 0;
        bool isDecomposed = false;

        size_t iIterations = 0;
        size_t iMaxIterations = -1;

        while (iabs(iDeterminant) ==0)
        {

            pMat = createMatrix(m_matN);
            LU = createMatrix(m_matN);

            pivot = createVec(m_matN);

            getRandomMatrix(pMat, m_matN);

            memcpy(LU, pMat, m_matN*m_matN*sizeof(int8_t));

            isDecomposed = decomposeMatrixRecipe(LU, pivot, m_matN);
            //int8_t d = 0;

            if (isDecomposed)
            {

                //applyRing(LU, n*n);

                iDeterminant = this->triangularDeterminant(LU,m_matN) % 2;// * determinant(LU,n);

                if (iDeterminant == 0)
                {
                    free(pMat);
                    free(LU);
                    free(pivot);
                }

            }

            ++iIterations;

        }

        if (iIterations >= iMaxIterations)
        {
            throw BijectiveMapping1DeterminantException(this);
        }

        // we have found our matrix, now invert

        int8_t* L = createMatrix(m_matN);
        int8_t* U = createMatrix(m_matN);

        for(size_t i=0; i < m_matN; ++i)
            for (size_t j = i; j < m_matN; ++j)
                U[pos(i,j)] = LU[pos(i,j)];

        for(size_t i=0; i < m_matN; ++i)
        {
            L[pos(i,i)] = 1;
            for (size_t j = 0; j < i; ++j)
                L[pos(i,j)] = LU[pos(i,j)];
        }

        free(pMat);

        int8_t* plumat = matrixMul(L, U, n);
        free(L);
        free(U);

        int8_t* nopiv = createVec(n);
        for (size_t i = 0; i < n; ++i)
            nopiv[i] = (int8_t) i;

        int8_t* inv = invertComb(LU, nopiv, n);
        free(LU);
        free(nopiv);

        applyRing(inv, n, n);

        int8_t* prod = matrixMul(plumat, inv, n);
        applyRing(prod, n, n);

        bool isValid = true;
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {

                if ( i == j)
                {

                    if (!(prod[pos(i,j)] == 1))
                    {
                        isValid = false;
                    }

                } else {

                    if (!(prod[pos(i,j)] == 0))
                    {
                        isValid = false;
                    }
                }

            }
        }

        if (isValid)
        {
            memcpy(pA, plumat, sizeof(int8_t) * n * n);
            memcpy(pAi, inv, sizeof(int8_t) * n * n);
        }

        free(plumat);
        free(inv);

        return isValid;

    }


    inline size_t npos(int i, int j, int n)
    {
        return i * n +j;
    }

    inline size_t pos(int i, int j)
    {
        return i * m_matN +j;
    }



    template <class T>
    void printVec(T* pVec)
    {
        std::cout << "[ ";

        for (uint32_t i =0; i < m_matN; ++i)
        {
            std::cout << std::to_string(pVec[i]) << " ";
        }

        std::cout << "]" << std::endl;
    }

    int8_t* createMatrix(size_t n)
    {
        return (int8_t*) calloc( sizeof(int8_t), n*n);
    }


    const uint32_t m_iK;
    const uint32_t m_matN;

    MemoryPool<FIELDTYPE>* m_pPool;

    int8_t* m_pA;
    int8_t* m_pAi;

    UBigInt* m_pArows;
    UBigInt* m_pAirows;

};


#endif //TSXCOUNT_BIJECTIVEKMAPPING_H
