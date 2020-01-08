//
// Created by mjopp on 05/12/2018.
//

#ifndef TSXCOUNT_MEMORYPOOL_H
#define TSXCOUNT_MEMORYPOOL_H

#include <stdlib.h>
#include <inttypes.h>
#include <omp.h>
#include <vector>
#include <iostream>
#include <thread>

template <class T>
class MemoryPool;

template <class T>
class MemoryPoolException: public std::exception
{
public:
    MemoryPoolException(MemoryPool<T>* pMapping, std::string sErrorMessage)
            : std::exception()
    {
        m_pMapping = pMapping;
        m_sErrorMessage = sErrorMessage;
    }


    virtual const char* what() const throw()
    {
        return this->getExceptText();
    }

    virtual const char* getExceptText() const
    {
        std::stringstream oSS;

        oSS << m_sErrorMessage;

        return oSS.str().c_str();
    }

protected:

    MemoryPool<T>* m_pMapping = NULL;
    std::string m_sErrorMessage;

};

struct MemLoc
{
    bool free;
    size_t ilength;
    void* addr;
};

template <class T>
class MemoryPatch
{
public:

    MemoryPatch(size_t elements)
    {
        this->m_pHeap = (T*) ::malloc(sizeof(T)* elements);
        MemLoc elem;
        elem.free = true;
        elem.ilength = elements;
        elem.addr = m_pHeap;

        // to avoid resize
        this->m_vHeap.reserve(elements);
        this->m_vHeap.push_back(elem);
        this->m_iInitialSize = elements;
    }

    ~MemoryPatch()
    {
        ::free(this->m_pHeap);
    }

    virtual MemLoc malloc(size_t iElemCount)
    {
        for (std::vector<MemLoc>::iterator oit = m_vHeap.begin(); oit != m_vHeap.end(); ++oit)
        {
            if ((oit->free) && (oit->ilength >= iElemCount))
            {

                MemLoc newelem;
                newelem.free = false;
                newelem.ilength = iElemCount;
                newelem.addr = oit->addr;

                // update old elem
                oit->addr = ((T*)oit->addr) + iElemCount;
                oit->ilength = oit->ilength - iElemCount;

                if (oit->ilength == 0)
                {
                    this->eraseIt(oit);
                }

                std::vector<MemLoc>::iterator oins = m_vHeap.insert(oit, newelem);

                //std::cout << m_vHeap[m_vHeap.size()-1].ilength << " " << m_vHeap[m_vHeap.size()-1].addr << std::endl;

                //long d = std::distance(m_vHeap.begin(), oins);

                return newelem;
            }
        }


        size_t iCnt=0;
        size_t iFreeCount = 0;

        for (auto elem : this->m_vHeap)
        {
            iCnt += elem.ilength;

            if (elem.free)
            {
                iFreeCount += elem.ilength;
            }
        }

        std::cerr << omp_get_thread_num() << " " << std::this_thread::get_id() << " Total element count " << iCnt << " with free spots: " << iFreeCount << std::endl;

        throw new MemoryPoolException<T>(NULL, "Unable to allocate memory.");

    }

    void eraseIt(std::vector<MemLoc>::iterator& oit)
    {
        long d = std::distance(m_vHeap.begin(), oit);
        m_vHeap.erase(oit);
    }

    void printLoadOut()
    {
        size_t iCnt=0;
        for (auto elem : this->m_vHeap)
        {
            iCnt += elem.ilength;
        }

        if (iCnt != this->m_iInitialSize)
        {
            std::cout << "error" << std::endl;
            std::cerr << omp_get_thread_num() << " " << std::this_thread::get_id() << " LOADOUT Total element count " << iCnt << std::endl;
        }

    }

    virtual void free(MemLoc& oElem)
    {

        for (std::vector<MemLoc>::iterator oit = m_vHeap.begin(); oit != m_vHeap.end(); ++oit)
        {
            if ((oit->free == false) && (oit->addr == oElem.addr))
            {
                long d = std::distance(m_vHeap.begin(), oit);

                if ((d+1 < m_vHeap.size()) && (m_vHeap[d+1].free == true))
                {

                    //std::cerr << "Freeing memory " << oElem.addr << std::endl;

                    m_vHeap[d+1].addr = oElem.addr;
                    m_vHeap[d+1].ilength += oElem.ilength;

                    this->eraseIt(oit);
                    this->printLoadOut();

                    return;
                } else if ((d >= 1) && (m_vHeap[d-1].free == true)) {

                    //std::cerr << "Freeing memory " << oElem.addr << std::endl;
                    m_vHeap[d-1].ilength += oElem.ilength;
                    this->eraseIt(oit);
                    this->printLoadOut();

                    return;
                } else {

                    //std::cerr << "Freeing memory " << oElem.addr << std::endl;
                    m_vHeap[d].free = true;
                    this->printLoadOut();
                    return;
                }
            }
        }

        throw new MemoryPoolException<T>(NULL, "Could not free Memory");


    }

protected:

    T* m_pHeap;
    std::vector<MemLoc> m_vHeap;
    size_t m_iInitialSize;
};

template <class T>
class NewAllocator : public MemoryPatch<T>
{
public:
    NewAllocator()
    : MemoryPatch<T>(0)
    {

    }

    virtual MemLoc* malloc(size_t iElemCount)
    {

        MemLoc oLoc;
        oLoc.free = true;
        oLoc.ilength = iElemCount;
        oLoc.addr = ::malloc(iElemCount * sizeof(T));

        this->m_lHeap.push_back(oLoc);

        MemLoc& lElem = this->m_lHeap.back();

        return &lElem;
    }

    virtual void free(MemLoc* pElem)
    {
        ::free(pElem);
    }
};


template <class T>
class MemoryPool {

public:
    MemoryPool(size_t iThreads)
    {
        m_iThreads = iThreads;
        m_ppPatches = (MemoryPatch<T>**) ::malloc(sizeof(MemoryPatch<T>*) * m_iThreads);

        for (size_t i = 0; i < m_iThreads; ++i)
        {
            m_ppPatches[i] = new MemoryPatch<T>(m_iDefaultSize);
        }
    }

    ~MemoryPool()
    {
        ::free(m_ppPatches);
    }

    void setThreads(size_t iThreads)
    {

        if (iThreads < m_iThreads)
        {
            return;
        }

        MemoryPatch<T>** oldPatches = m_ppPatches;

        m_ppPatches = (MemoryPatch<T>**) ::malloc(sizeof(MemoryPatch<T>*) * iThreads);

        size_t i;
        for (i = 0; i < m_iThreads; ++i)
        {
            m_ppPatches[i] = oldPatches[i];
        }

        for (; i < iThreads; ++i)
        {
            m_ppPatches[i] = new MemoryPatch<T>(m_iDefaultSize);
        }

        m_iThreads = iThreads;

    }

    virtual MemLoc malloc(size_t iElemCount)
    {
        size_t iThreadNum = omp_get_thread_num();
        return m_ppPatches[iThreadNum]->malloc(iElemCount);
    }

    virtual void free(MemLoc& pElem)
    {
        size_t iThreadNum = omp_get_thread_num();
        m_ppPatches[iThreadNum]->free(pElem);
    }

protected:

    size_t m_iThreads;
    MemoryPatch<T>** m_ppPatches;

    const size_t m_iDefaultSize = 200000;

};


template <class T>
class DirectMemoryPool : public MemoryPool<T>{

public:
    DirectMemoryPool()
    : MemoryPool<T>(0)
    {

    }

    ~DirectMemoryPool()
    {
    }

    MemLoc malloc(size_t iElemCount)
    {

        MemLoc oRet;
        oRet.ilength = iElemCount;
        oRet.free = false;
        oRet.addr = ::malloc(sizeof(T) * iElemCount);

        return oRet;
    }

    void free(MemLoc& pElem)
    {
        ::free(pElem.addr);
    }

protected:


};




#endif //TSXCOUNT_MEMORYPOOL_H
