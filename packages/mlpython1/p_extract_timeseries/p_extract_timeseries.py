import sys
import numpy as np
import os
import traceback

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# imports from mlpy
from mlpy import ProcessorManager
from mlpy import DiskReadMda, readmda, writemda32, DiskWriteMda
from timeserieschunkreader import TimeseriesChunkReader

class extract_timeseries:
    """
    Extract a chunk of a timeseries dataset and possibly a subset of channels

    Parameters
    ----------
    timeseries : INPUT
        Path of timeseries, MxN where M is number of channels and N number of timepoints, in either .mda or raw binary format. If raw binary, then you must supply dtype and num_channels.
    channels_array : INPUT 
        (optional) Path of array of channel numbers (positive integers). Either use this or the channels parameter, not both.
        
    timeseries_out : OUTPUT
        Path of output timeseries in .mda format    
        
    channels : string
        (Optional) Comma-separated list of channels to extract. Either use this or the channels_array input, not both.
    t1 : integer
        (Optional) Integer start timepoint (zero-based indexing). If -1 will set to zero.
    t2 : integer
        (Optional) Integer end timepoint (zero-based indexing). If -1 will set to N-1."},
    timeseries_dtype : string
        (Optional) Only supply this if timeseries is in raw binary format. Choices are int16, uint16, int32, float32, etc.
    timeseries_num_channels : integer
        (Optional) Only supply this if timeseries is in raw binary format. Integer representing number of channels. Number of timepoints will be deduced
    """    
    name='mlpython1.extract_timeseries'
    version="0.1"
    def __call__(self,*,
            timeseries,
            channels_array='',
            timeseries_out,
            channels='',t1=-1,t2=-1,
            timeseries_dtype='',timeseries_num_channels=0
        ):
        if channels:
            self._channels=np.fromstring(channels,dtype=int,sep=',')
        elif channels_array:
            self._channels=channels_array
        else:
            self._channels=np.empty(0)
            
        t1=int(t1)
        t2=int(t2)
        
        X=DiskReadMda(timeseries)
        M,N = X.N1(),X.N2()
        if (self._channels.size==0):
            self._channels=np.array(1+np.arange(M))
        M2=self._channels.size
        
        if (t1<0):
            t1=0
        if (t2<0):
            t2=N-1
            
        self._writer=DiskWriteMda(timeseries_out,[M2,N])

        num_mb=10
        chunk_size=np.maximum(100,int(num_mb*1e6/(M*4)))
        chunk_size=10
        TCR=TimeseriesChunkReader(chunk_size=chunk_size, overlap_size=0, t1=t1, t2=t2)
        return TCR.run(timeseries,self._kernel)            
    def _kernel(self,chunk,info):
        chunk=chunk[(self._channels-1).tolist(),]
        return self._writer.writeChunk(chunk,i1=0,i2=info.t1)
    def test(self):
        try:
            M,N = 4,100
            X=np.random.rand(M,N)
            writemda32(X,'tmp.mda')
            ret=self(timeseries="tmp.mda",timeseries_out="tmp2.mda",channels="1,3",t1="-1",t2="-1")
            assert(ret)
            A=readmda('tmp.mda')
            B=readmda('tmp2.mda')
            assert(B.shape[0]==2)
            assert(B.shape[1]==N)
            assert(np.array_equal(A[[0,2],],B))
            return True 
        except Exception as e:
            traceback.print_exc()
            return False

#import numpydoc
#doc=numpydoc.docscrape.FunctionDoc(extract_timeseries2)
#params=doc["Parameters"]
#for j in range(len(params)):
#    print(params[j])
#    print("")

#P=extract_timeseries()
#ret=P.test()
#print ("Test result: %d" % (ret))

PM=ProcessorManager()
PM.registerProcessor(extract_timeseries)
if not PM.run(sys.argv):
    exit(-1)