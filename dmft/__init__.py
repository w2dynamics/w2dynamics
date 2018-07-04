"""
Package providing DMFT stuff

MPI and modern programming
--------------------------
(A somewhat opiniated essay describing the design criteria for the MPI parts)

MPI, the powerful, low-level de-facto standard for massive parallelisation,
strongly encourages the use of the master/slave pattern: there is a "master" or
root process, which distributes work (MPI `scatter`, `bcast`) and data to the
possibly interwoven "slave" processes and collects their results (MPI `gather`,
`reduce`) at the end. This is fine for unstructured code, however it begins to
feel odd as soon as one moves to procedural code and breaks down for functional
codes.

         Pre-function broadcast                Inner-function broadcast
                                              __________________________
         ____    _____________               |   ____    _____________  |
    x ->| r0 |->|             |-> f      x ---->| r0 |->|             |---> f
        |    |  | MPI-enabled |              |  |    |  | MPI-enabled | |
        | r1 |->|     f(x)    |-> ?      ? ->|? | r1 |->| computation |---> ?
        |    |  |             |              |  |    |  |             | |
        :    :  :             :              |  :    :  :             : :
        | rn |->|             |-> ?      ? ->|? | rn |->|             |---> ?
        |____|  |_____________|              |  |____|  |_____________| |
                                             |__________________________|

Consider an arbitrary function that performs a MPI-parallelised calculation:
every process needs to call the function in order to be able to run it in
parallel. In the master/slave model, this means that the data has to be
distributed either before the function call or inside the function. In both
cases we end up with different function signatures for the master and slave
processes, because either the return code or the parameter list are meaningless
on the slaves. Also, pre-function distribution is a user pitfall and makes the
caller code unnecessarily MPI-aware, while with in-function distribution, it
seems unclear on which nodes parameters and return codes are filled with
meaningful values.

While some flaws can be remedied with comprehensive interface documentation,
they inherently break functional code. Suppose we want to supply our function
as argument to a third-party ODE driver routine: now *both* approaches fail,
since we cannot modify the caller code and cannot commu-nicate on which ranks
the result or the parameters are meaningful, so the ODE routine might call the
routine on different cores with different parameters, or even a different
number of times. A solution -- apart from writing your own solver, which
sadly happens often enough -- is to abandon the master/slave concept in favour
of a "consistent" interface: the function expects the identical parameters and
returns identical values for every process (cf. MPI `alltoall`, `allgather`,
`allreduce`):


        Broadcast execution                 MPI consistent interface model
         ________________                          ________________
    x ->| f(x) on rank 0 |-> f                x ->|                |-> f
        |________________|                        |  MPI-enabled   |
    x ->| f(x) on rank 1 |-> f                x ->|     f(x)       |-> f
        |________________|                        |                |
         _______:________                         :                :
    x ->| f(x) on rank n |-> f                x ->|                |-> f
        |________________|                        |________________|


The consistent interface allows to hide MPI as implementation detail inside the
function, which means MPI-enabled functions can be used as drop-in replacements
for single-core defaults. Furthermore, it removes the necessity of master/slave
branching instructions, which makes the programme "look" single-core. So why is
it not more widely used? Aside from being lesser known, there are three main
drawbacks:

  1. there is a performance overhead, because results have to be distributed to
     every core and code execution between functions is repeated ("broadcast")
     on any processes;

  2. it fails when the memory demand exceeds the resources of one core; and

  3. it duplicates (sequential) I/O over the cores, which is a bottleneck for
     read operations and yields duplication for write operations.

These problems can in principle be solved by moving to a dispatcher model,
where the main branch code is executed only on the master process, which then
distributes work packages to the slave processes. Dispatchers are difficult to
implement and to generalise, hard to optimize and not necessarily the ideal MPI
topology, so we will try to design around these flaws:

Typically, performance is not an issue unless the MPI function is used in a
very tight loop or a significant runtime share is spent in process-broadcast
execution. For the memory issue, one can use an encapsulation pattern: the
memory is distributed over several processes and a function is used to
retrieve, set, or alter parts of it. One need not provide a general-purpose
interface, but only implement the use cases in the program. Sequential I/O
operations need to be guarded against parallel execution. Since in scientific
applications, such I/O is typically performed by the main programme, it is fine
to use master/slave branching there.
"""
