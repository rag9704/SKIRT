#include "ParallelTable.hpp"
#include "FatalError.hpp"
#include "PeerToPeerCommunicator.hpp"
#include "ProcessAssigner.hpp"
#include "TimeLogger.hpp"

ParallelTable::ParallelTable()
    : _name(""), _colAssigner(0), _rowAssigner(0), _distributed(false), _synced(false), _initialized(false), _comm(0), _log(0)
{
}

////////////////////////////////////////////////////////////////////

void ParallelTable::initialize(QString name, const ProcessAssigner *colAssigner, const ProcessAssigner *rowAssigner, writeState writeOn)
{
    _name = name;
    _colAssigner = colAssigner;
    _rowAssigner = rowAssigner;
    _writeOn = writeOn;

    _comm = colAssigner->find<PeerToPeerCommunicator>();
    _log = colAssigner->find<Log>();

    _totalRows = _rowAssigner->total();
    _totalCols = _colAssigner->total();

    // Use the distributed memory scheme
    if (colAssigner->parallel() && rowAssigner->parallel() && _comm->isMultiProc())
    {
        _log->info(_name + " is distributed. Sizes of local tables are ("
                   + QString::number(_totalRows) + "," + QString::number(_colAssigner->nvalues())
                   + ") and ("
                   + QString::number(_rowAssigner->nvalues()) + "," + QString::number(_totalCols) + ")"
                   );

        _distributed = true;
        _columns.resize(_totalRows,colAssigner->nvalues());
        _rows.resize(rowAssigner->nvalues(),_totalCols);

        int Nprocs = _comm->size();
        _displacementvv.reserve(Nprocs);
        for (int r=0; r<Nprocs; r++) _displacementvv.push_back(_colAssigner->indicesForRank(r));
    }
    // not distributed
    else
    {
        _log->info(_name + " is not distributed. Size is ("
                    + QString::number(_totalRows) + "," + QString::number(_totalCols) + ")"
                   );

        _distributed = false;
        if (writeOn == COLUMN) _columns.resize(_totalRows,_totalCols);
        else if (writeOn == ROW) _rows.resize(_totalRows,_totalCols);
        else throw FATALERROR("Invalid writeState for ParallelTable");
    }
    _initialized = true;
    _synced = true;
}

////////////////////////////////////////////////////////////////////

const double& ParallelTable::operator()(size_t i, size_t j) const
{
    if (!_synced) throw FATALERROR(_name + " says: sync() must be called before using the read operator");

    // WORKING DISTRIBUTED: Read from the table opposite to the table we _writeOn.
    if (_distributed)
    {
        if (_writeOn == COLUMN) // read from _rows
        {
            if (!_rowAssigner->validIndex(i))
                throw FATALERROR(_name + " says: Row of ParallelTable not available on this process");
            return _rows(_rowAssigner->relativeIndex(i),j);
        }
        else                    // read from _columns
        {
            if (!_colAssigner->validIndex(j))
                throw FATALERROR(_name + " says: Column of ParallelTable not available on this process");
            return _columns(i,_colAssigner->relativeIndex(j));
        }
    }
    // WORKING NON-DISTRIBUTED: Reading and writing happens in the same table.
    else
    {
        if (_writeOn == COLUMN)
            return _columns(i,j);
        else
            return _rows(i,j);
    }
}

////////////////////////////////////////////////////////////////////

double& ParallelTable::operator()(size_t i, size_t j)
{
    _synced = false;

    // WORKING DISTRIBUTED: Writable reference to the table we _writeOn.
    if (_distributed)
    {
        if (_writeOn == COLUMN) // Write on _columns
        {
            if (!_colAssigner->validIndex(j))
                throw FATALERROR(_name + " says: Column of ParallelTable not available on this process");
            return _columns(i,_colAssigner->relativeIndex(j));
        }
        else                    // Write on _rows
        {
            if (!_rowAssigner->validIndex(i))
                throw FATALERROR(_name + " says: Row of ParallelTable not available on this process");
            return _rows(_rowAssigner->relativeIndex(i),j);
        }
    }
    // WORKING NON-DISTRIBUTED: Reading and writing happens in the same table.
    else
    {
        if (_writeOn == COLUMN)
            return _columns(i,j);
        else
            return _rows(i,j);
    }
}

////////////////////////////////////////////////////////////////////

Array& ParallelTable::operator[](size_t i)
{
    _synced = false;

    if (_writeOn == ROW) // return writable reference to a complete row
    {
        if (!_rowAssigner->validIndex(i))
            throw FATALERROR(_name + " says: Row of ParallelTable not available on this process");
        else return _rows[_distributed ? _rowAssigner->relativeIndex(i) : i];
    }
    else throw FATALERROR(_name + " says: This ParallelTable does not have writable rows.");
}

////////////////////////////////////////////////////////////////////

const Array& ParallelTable::operator[](size_t i) const
{
    if (!_synced) throw FATALERROR(_name + " says: sync() must be called before asking a read reference");

    if (_distributed && _writeOn == COLUMN) // return read only reference to a complete row
    {
        if (!_rowAssigner->validIndex(i))
            throw FATALERROR(_name + " says: Row of ParallelTable not available on this process");
        else return _rows[_rowAssigner->relativeIndex(i)];
    }
    else throw FATALERROR(_name + " says: This ParallelTable does not have readable rows.");
}

////////////////////////////////////////////////////////////////////

void ParallelTable::sync()
{
    if (!_synced)
    {
        TimeLogger logger(_log->verbose() && _comm->isMultiProc() ? _log : 0, "communication of " + _name);

        if (!_distributed) sum_all();
        else if (_writeOn == COLUMN) experimental_col_to_row();
        else if (_writeOn == ROW) experimental_row_to_col();
    }
    _synced = true;
}

////////////////////////////////////////////////////////////////////

void ParallelTable::clear()
{
    _columns.clear();
    for (size_t i=0; i<_rows.size(0); i++) _rows[i] *= 0;

    _synced = true;
}

////////////////////////////////////////////////////////////////////

double ParallelTable::sumRow(size_t i) const
{
    if (!_synced) throw FATALERROR(_name + " says: sync() must be called before using summation functions");

    double sum = 0;
    if (_distributed)
    {
        if (_rowAssigner->validIndex(i)) // we have the whole row
            return _rows[_rowAssigner->relativeIndex(i)].sum();
        else throw FATALERROR(_name + " says: sumRow(index) called on wrong index");
    }
    else // not distributed -> everything available, so straightforward
    {
        for (int j=0; j<_totalCols; j++)
            sum += (*this)(i,j);
        return sum;
    }
}

////////////////////////////////////////////////////////////////////

double ParallelTable::sumColumn(size_t j) const
{
    if (!_synced) throw FATALERROR(_name + " says: sync() must be called before using summation functions");

    double sum = 0;
    if (_distributed)
    {
        if(_colAssigner->validIndex(j)) // we have the whole column
            for (int i=0; i<_totalRows; i++)
                sum += _columns(i,_colAssigner->relativeIndex(j));
        else throw FATALERROR(_name + " says: sumColumn(index) called on wrong index");
    }
    else // not distributed
    {
        for (int i=0; i<_totalRows; i++)
            sum += (*this)(i,j);
    }
    return sum;
}

////////////////////////////////////////////////////////////////////

Array ParallelTable::stackColumns() const
{
    if (!_synced) throw FATALERROR(_name + " says: sync() must be called before using summation functions");

    Array result(_totalRows);

    if(_distributed)
    {
        // fastest way: sum only the values in _rows, then do one big sum over the processes
        for (size_t iRel=0; iRel<_rows.size(0); iRel++)
            result[_rowAssigner->absoluteIndex(iRel)] = _rows[iRel].sum();

        TimeLogger logger(_log->verbose() ? _log : 0, "summing columns in " + _name);

        _comm->sum_all(result);
    }
    else
    {
        for (size_t i=0; i<_rowAssigner->total(); i++)
            result[i] = sumRow(i);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

Array ParallelTable::stackRows() const
{
    if (!_synced) throw FATALERROR(_name + " says: sync() must be called before using summation functions");

    Array result(_totalCols);

    if(_distributed)
    {
        for (size_t jRel=0; jRel<_columns.size(1); jRel++)
        {
            size_t j = _colAssigner->absoluteIndex(jRel);
            for (size_t i=0; i<_columns.size(0); i++)
                result[j] += _columns(i,jRel);
        }

        TimeLogger logger(_log->verbose() ? _log : 0, "summing rows in " + _name);

        _comm->sum_all(result);
    }
    else
    {
        for (int j=0; j<_totalCols; j++)
            result[j] = sumColumn(j);
    }
    return result;
}

////////////////////////////////////////////////////////////////////

double ParallelTable::sumEverything() const
{
    if (!_synced) throw FATALERROR(_name + " says: sync() must be called before using summation functions");

    return stackColumns().sum();
}

////////////////////////////////////////////////////////////////////

bool ParallelTable::distributed() const
{
    return _distributed;
}

////////////////////////////////////////////////////////////////////

bool ParallelTable::initialized() const
{
    return _initialized;
}

////////////////////////////////////////////////////////////////////

void ParallelTable::sum_all()
{
    if (_writeOn == COLUMN)
    {
        Array& arr = _columns.getArray();
        _comm->sum_all(arr);
    }
    else
    {
        for (int i=0; i<_totalRows; i++)
        {
            Array& arr = _rows[i];
            _comm->sum_all(arr);
        }
    }
}

////////////////////////////////////////////////////////////////////

void ParallelTable::col_to_row()
{
    int thisRank = _comm->rank();

    int sendCount = 0;
    int N = 750;

    for (int i=0; i<_totalRows; i++) // for each possible row of the big array (determines receiver)
    {
        int tgtRank = _rowAssigner->rankForIndex(i); // the rank where target has row i in it

        // post some receives
        for (int j=0; j<_totalCols; j++) // for each possible column of the big array (determines sender)
        {
            int srcRank = _colAssigner->rankForIndex(j); // the rank where the source has col j in it
            int tag = i*_totalCols + j; // unique for each position in the big array

            if (thisRank!= srcRank && thisRank == tgtRank) // receive if target needs row i
            {
                double& recvbuf = _rows(_rowAssigner->relativeIndex(i),j);
                _comm->receiveDouble(recvbuf,srcRank,tag);
            }
        }
        // post some sends
        for (int j=0; j<_totalCols; j++) // for each possible column of the big array (determines sender)
        {
            int srcRank = _colAssigner->rankForIndex(j); // the rank where the source has col j in it
            int tag = i*_totalCols + j; // unique for each position in the big array

            if (thisRank == srcRank) // if this process has column j, do send
            {
                double& sendbuf = _columns(i,_colAssigner->relativeIndex(j));
                if (thisRank == tgtRank) _rows(_rowAssigner->relativeIndex(i), j) = sendbuf;
                else _comm->sendDouble(sendbuf,tgtRank,tag);
            }
            sendCount++;
        }
        // clear the request vector every N
        if (sendCount > N)
        {
            _comm->finishRequests();
            sendCount = 0;
        }
    }
    _comm->finishRequests();
}

////////////////////////////////////////////////////////////////////

void ParallelTable::row_to_col()
{
    int thisRank = _comm->rank();

    int sendCount = 0;
    int N = 750;

    for (int i=0; i<_totalRows; i++) // for each row (determines sender)
    {
        int srcRank = _rowAssigner->rankForIndex(i);

        // post some receives
        for (int j=0; j<_totalCols; j++) // for each column (determines receiver)
        {
            int tgtRank = _colAssigner->rankForIndex(j);
            int tag = i*_totalCols + j;

            if (thisRank != srcRank && thisRank == tgtRank)
            {
                double& recvbuf = _columns(i,_colAssigner->relativeIndex(j));
                _comm->receiveDouble(recvbuf,srcRank,tag);
            }
        }
        // post some sends
        for (int j=0; j<_totalCols; j++) // for each column (determines receiver)
        {
            int tgtRank = _colAssigner->rankForIndex(j);
            int tag = i*_totalCols + j;

            if (thisRank == srcRank)
            {
                double& sendbuf = _rows(_rowAssigner->relativeIndex(i),j);
                if (thisRank == tgtRank) _columns(i,_colAssigner->relativeIndex(j)) = sendbuf;
                else _comm->sendDouble(sendbuf,tgtRank,tag);
            }
            sendCount++;
        }
        // clear the request vector every N
        if (sendCount > N)
        {
            _comm->finishRequests();
            sendCount = 0;
        }
    }
    _comm->finishRequests();
}

////////////////////////////////////////////////////////////////////

void ParallelTable::experimental_col_to_row()
{
    // prepare to receive using doubles according to the given displacements for each process
    _comm->presetConfigure(1,_displacementvv);

    // for each row
    for (int i=0; i<_totalRows; i++)
    {
        double* sendBuffer = &_columns(i,0);
        double* recvBuffer = _rowAssigner->validIndex(i) ? &_rows(_rowAssigner->relativeIndex(i),0) : 0;
        int recvRank = _rowAssigner->rankForIndex(i);

        _comm->presetGatherw(sendBuffer, _columns.size(1), recvBuffer, recvRank);

        //_comm->gatherw(sendBuffer, _columns.size(1), recvBuffer, recvRank, 1, _displacementvv);
        // Elk proces i zendt vanuit sendBuffer, _columns.size(1) doubles naar recvBuffer op recvRank in de vorm
        // van datatype i bepaald door blocksize 1 en _displacementvv[i].
    }

    // free the memory taken by the datatypes
    _comm->presetClear();
}

////////////////////////////////////////////////////////////////////

void ParallelTable::experimental_row_to_col()
{
    // prepare to send using doubles according to the given displacements for each process
    _comm->presetConfigure(1,_displacementvv);

    for (int i=0; i<_totalRows; i++)
    {
        double* sendBuffer = _rowAssigner->validIndex(i) ? &_rows(_rowAssigner->relativeIndex(i),0) : 0;
        double* recvBuffer = &_columns(i,0);
        int sendRank = _rowAssigner->rankForIndex(i);

        _comm->presetScatterw(sendBuffer, sendRank, recvBuffer, _columns.size(1));

        //_comm->scatterw(sendBuffer, sendRank, 1, _displacementvv, recvBuffer, _columns.size(1));
        // Proces sendRank zendt vanuit sendBuffer een aantal doubles in de vorm van datatype i bepaald door
        // _displacementvv[i] naar recvBuffer op rank i.
    }

    _comm->presetClear();
}

////////////////////////////////////////////////////////////////////