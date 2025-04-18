# from mpi4py import MPI
# import numpy as np

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()
# N = 40
# print(f"Hellp from rank {rank} / {size}")
# print("rank = ",rank)
# A = np.zeros(N,dtype='d')
# for i in range(rank, N, size):
#     print(i)
#     A[i] = i
# total = None
# if rank==0:
#     total = np.empty((N,),dtype='d')
# comm.Gather(A,total,root=0)
# print("A = ",A)
# print("total = ",total)
from mpi4py import MPI
import numpy

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# passing MPI datatypes explicitly
if rank == 0:
    data = numpy.arange(1000, dtype='i')
    comm.Send([data, MPI.INT], dest=1, tag=77)
elif rank == 1:
    data = numpy.empty(1000, dtype='i')
    comm.Recv([data, MPI.INT], source=0, tag=77)

# automatic MPI datatype discovery
if rank == 0:
    data = numpy.arange(10, dtype=numpy.float64)
    comm.Send(data, dest=1, tag=13)
elif rank == 1:
    data = numpy.empty(10, dtype=numpy.float64)
    comm.Recv(data, source=0, tag=13)

print(data)