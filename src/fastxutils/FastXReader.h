//
// Created by mjoppich on 11/1/17.
//

#ifndef TSXCOUNT_FASTXREADER_H
#define TSXCOUNT_FASTXREADER_H

#include <string>
#include <iostream>
#include <zlib.h>
#include <vector>
#include <fstream>
#include <utils/CLParser.h>
#include <utils/Utils.h>

#define FASTQREADER_ZLIB_CHUNK 65536

class FASTXreaderException: public std::exception
{
public:
    FASTXreaderException(std::string sText)
            : std::exception(), m_sText(sText)
    {
    }


    virtual const char* what() const throw()
    {
        return m_sText.c_str();
    }


protected:

    const std::string m_sText;

};

class IParseEntry
{
public:

    std::string getID()
    {
        return sReadID;
    }

    std::string getSequence()
    {
        return sSequence;
    }

    virtual int getLinesRequired() const = 0;

protected:

    std::string sReadID;
    std::string sSequence;

};

class FASTQEntry : public IParseEntry {
public:

    FASTQEntry()
    :IParseEntry()
    {

    }

    FASTQEntry(std::vector<std::string>& lines)
    {
        sReadID   = lines[0];
        sSequence = lines[1];
        // +
        sQuality  = lines[3];

        /*
        if (sReadID[0] != '@')
        {
            std::cout << "bla" << std::endl;
        }
         */
    }

    int getLinesRequired() const
    {
        return 4;
    }


protected:
    std::string sQuality;

};

class FASTAEntry : public IParseEntry {
public:

    FASTAEntry()
            :IParseEntry()
    {

    }

    FASTAEntry(std::vector<std::string>& lines)
    {
        sReadID   = lines[0];
        sSequence = lines[1];
    }

    int getLinesRequired() const
    {
        return 2;
    }
};

template <class T >
class FASTXreader {

private:
    void ValidateType( FASTQEntry    &i ) const {}
    void ValidateType( FASTAEntry    &i ) const {}


public:

    static FASTXreader<IParseEntry>* createReader(CLParser* pConfig)
    {

        if (pConfig->isSet("fastq") || pConfig->isSet("fq"))
        {
            return (FASTXreader<IParseEntry>*) FASTXreader::createFQReader(pConfig);
        }

        if (pConfig->isSet("fasta") || pConfig->isSet("fa"))
        {
            return (FASTXreader<IParseEntry>*) FASTXreader::createFQReader(pConfig);
        }

    }

    static FASTXreader<FASTQEntry>* createFQReader(CLParser* pConfig)
    {

        if (pConfig == NULL)
            return NULL;

        std::string* pFASTQFile = pConfig->getFileName("fastq");

        bool bError = false;
        if (pFASTQFile == NULL)
        {
            std::string* pFASTQFile = pConfig->getFileName("fq");

            if (pFASTQFile == NULL)
            {
                bError = true;
            }
        }

        if (bError)
        {

            pConfig->printAllArguments();
            throw new FASTXreaderException("No input file found for config");

            return NULL;

        }

        std::cout << "FASTQ: " << *pFASTQFile << std::endl;

        return new FASTXreader<FASTQEntry>(pFASTQFile);

    }

    FASTXreader(std::string* pFileName, bool bZLIBMode = false)
            : m_bZLIBMode(bZLIBMode), m_pFileName(pFileName)
    {

        T valid;
        ValidateType( valid );

        std::string sGZ = ".gz";

        if ( pFileName->rfind( sGZ ) == pFileName->length()-3 )
        {
            m_bZLIBMode = true;
        }

        if (m_bZLIBMode)
        {
            m_pGZIPFile = fopen (m_pFileName->c_str(), "rb");
        } else {

            m_oInputStream.open(m_pFileName->c_str());

            if (m_oInputStream.eof() == true)
                return;

        }



    }

    void close()
    {

        if (m_bZLIBMode)
        {
            this->finishGZIPFile();
        } else {
            this->finishFile();
        }

    }


    void readEntries(std::vector<T>* pEntries, size_t iChunkSize = 1)
    {

        T elem;
        uint8_t iLinesPerEntry = elem.getLinesRequired();

        size_t iLinesNeeded = iChunkSize * iLinesPerEntry;
        this->readLinesFromFile(iLinesNeeded);

        if (m_vBufferedLines.size() < iLinesNeeded)
        {
            std::cerr << "FASTQreader: less lines available than requested" << std::endl;
        }

        iLinesNeeded = std::min(iLinesNeeded, m_vBufferedLines.size());

        pEntries->reserve( iLinesNeeded/ iLinesPerEntry );

        for (size_t i = 0; i < iLinesNeeded; i += iLinesPerEntry)
        {

            std::vector<std::string> entryLines;

            for (size_t j = 0; j < iLinesPerEntry; ++j)
            {
                entryLines.push_back(m_vBufferedLines[i+j]);
            }

            T oEntry( entryLines );

            pEntries->push_back(oEntry);

        }

        m_vBufferedLines.erase(m_vBufferedLines.begin(), m_vBufferedLines.begin()+iLinesNeeded);

    }


    std::vector<T> readEntries(size_t iChunkSize = 1)
    {

        std::vector<T> vEntries;

        this->readEntries(&vEntries, iChunkSize);

        return vEntries;

    }

    std::vector<T>* getEntries(size_t iChunkSize = 1)
    {

        std::vector<T>* pEntries = new std::vector<T>();

        this->readEntries( pEntries, iChunkSize);

        return pEntries;

    }

    bool hasNext()
    {
        if (m_vBufferedLines.size() > 0)
            return true;


        if (m_bZLIBMode)
        {

            return m_pGZIPFile != NULL;


        } else {

            return m_oInputStream.is_open() && !m_oInputStream.eof();

        }


        return false;
    }


protected:

    void readLinesFromFile(size_t iLines)
    {

        if (iLines > m_vBufferedLines.size())
        {
            // we need to fetch new lines

            if (m_bZLIBMode)
            {
                char pInBuffer[ FASTQREADER_ZLIB_CHUNK ];

                while (m_pGZIPFile != NULL) {

                    uint32_t iBytesRead = fread (pInBuffer, sizeof (char), sizeof (pInBuffer), m_pGZIPFile);

                    // THIS FETCHES ALL LINES!
                    std::string sZLIBChunk = this->zlib_uncompress(pInBuffer, iBytesRead);

                    m_sBufferedText.append( sZLIBChunk );
                    std::vector<std::string> vChunkLines = Utils::split(m_sBufferedText, '\n');

                    bool bLastWasNewline = m_sBufferedText.at(m_sBufferedText.size()-1) == '\n';
                    size_t iCopyLines = vChunkLines.size()-1;

                    if (bLastWasNewline)
                    {
                        ++iCopyLines;
                        m_sBufferedText = "";
                    } else {
                        m_sBufferedText = vChunkLines[vChunkLines.size()-1];
                    }

                    m_vBufferedLines.insert(m_vBufferedLines.end(), vChunkLines.begin(), vChunkLines.begin() + iCopyLines);


                    // if we are finished, this closes the file
                    if (feof (m_pGZIPFile)) {

                        finishGZIPFile();

                        break;
                    }

                    //std::cout << bLastWasNewline << " " << m_vBufferedLines.size() << std::endl;

                    if (m_vBufferedLines.size() > iLines)
                        break;
                }

            } else {


                std::string sLine;

                m_vBufferedLines.reserve(iLines);

                while ((!m_oInputStream.eof()) && ( m_vBufferedLines.size() < iLines )) {

                    std::getline(m_oInputStream, sLine);

                    if (sLine.length() == 0)
                        continue;

                    m_vBufferedLines.push_back(sLine);

                }

                if (m_oInputStream.eof())
                {
                    this->finishFile();
                }

            }



        }

    }

    std::string zlib_uncompress(char* pInBuffer, uint32_t iBytes)
    {
        int ret;
        uint32_t have;
        char out[FASTQREADER_ZLIB_CHUNK];

        if (m_pStream == NULL)
        {
            m_pStream = new z_stream();

            m_pStream->zalloc = Z_NULL;
            m_pStream->zfree = Z_NULL;
            m_pStream->opaque = Z_NULL;
            m_pStream->avail_in = 0;
            m_pStream->next_in = Z_NULL;

            ret = inflateInit2(m_pStream, 32 + MAX_WBITS);
            if (ret != Z_OK)
            {
                throw FASTXreaderException("inflateInit2 returned code " + std::to_string(ret));
            }
        }



        std::string sUncompressData = "";

        m_pStream->avail_in = iBytes;
        m_pStream->next_in = (unsigned char*) pInBuffer;
        do {
            m_pStream->avail_out = FASTQREADER_ZLIB_CHUNK;
            m_pStream->next_out = (unsigned char*) out;
            ret = inflate(m_pStream, Z_NO_FLUSH);

            if (ret != Z_OK)
            {
                switch (ret) {
                    case Z_NEED_DICT:
                    case Z_DATA_ERROR:
                    case Z_MEM_ERROR:
                    default:
                        inflateEnd(m_pStream);
                        throw FASTXreaderException("inflate returned code " + std::to_string(ret));
                }
            }


            have = FASTQREADER_ZLIB_CHUNK - m_pStream->avail_out;
            sUncompressData.append((char*)out, have);
        } while (m_pStream->avail_out == 0);


        return sUncompressData;
    }

    void finishGZIPFile()
    {
        int ret = inflateEnd(m_pStream);
        if (ret != Z_OK) {
            throw FASTXreaderException("inflateEnd returned code " + std::to_string(ret));
        }

        m_pStream = NULL;

        fclose(m_pGZIPFile);
        m_pGZIPFile = NULL;
    }

    void finishFile()
    {

        if (m_oInputStream.is_open())
            m_oInputStream.close();

    }


    bool m_bZLIBMode = false;
    std::string* m_pFileName;

    std::string m_sBufferedText = "";
    std::vector<std::string> m_vBufferedLines;

    // ZLIB mode
    FILE* m_pGZIPFile;
    z_stream* m_pStream = NULL;


    // file mode
    std::ifstream m_oInputStream;

};


#endif //TSXCOUNT_FASTXREADER_H
