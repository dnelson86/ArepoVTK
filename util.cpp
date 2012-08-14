/*
 * util.cpp
 * dnelson
 */
 
#include "util.h"

#include <fcntl.h>
#include <unistd.h>

// Timer
Timer::Timer()
{
    time0 = elapsed = 0;
    running = 0;
}

double Timer::GetTime()
{
    gettimeofday( &timeofday, NULL );
    return timeofday.tv_sec + timeofday.tv_usec / 1000000.0;
}

void Timer::Start()
{
    running = 1;
    time0 = GetTime();
}

void Timer::Stop()
{
    running = 0;
    elapsed += GetTime() - time0;
}

void Timer::Reset()
{
    running = 0;
    elapsed = 0;
}

double Timer::Time()
{
    if (running) {
        Stop();
        Start();
    }
    return elapsed;
}

// memory

void *AllocAligned(size_t size) {
    return memalign(L1_CACHE_LINE_SIZE, size);
}

void FreeAligned(void *ptr) {
    if (!ptr) return;
    free(ptr);
}

// random number generation

#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

void RNG::Seed(uint32_t seed) const {
    mt[0]= seed & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
        (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* generates a random number on [0,1)-real-interval */
float RNG::RandomFloat() const
{
    float v = (RandomUInt() & 0xffffff) / float(1 << 24);
    return v;
}

unsigned long RNG::RandomUInt() const
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if Seed() has not been called, */
            Seed(5489UL); /* default initial seed */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

void StratifiedSample1D(float *samp, int nSamples, RNG &rng, bool jitter) {
    float invTot = 1.f / nSamples;
    for (int i = 0;  i < nSamples; ++i) {
        float delta = jitter ? rng.RandomFloat() : 0.5f;
        *samp++ = min((i + delta) * invTot, OneMinusEpsilon);
    }
}

void StratifiedSample2D(float *samp, int nx, int ny, RNG &rng, bool jitter) {
    float dx = 1.f / nx, dy = 1.f / ny;
    for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x) {
            float jx = jitter ? rng.RandomFloat() : 0.5f;
            float jy = jitter ? rng.RandomFloat() : 0.5f;
            *samp++ = min((x + jx) * dx, OneMinusEpsilon);
            *samp++ = min((y + jy) * dy, OneMinusEpsilon);
        }
}

void LatinHypercube(float *samples, uint32_t nSamples, uint32_t nDim, RNG &rng) {
    // Generate LHS samples along diagonal
    float delta = 1.0f / nSamples;
    for (uint32_t i = 0; i < nSamples; ++i)
        for (uint32_t j = 0; j < nDim; ++j)
            samples[nDim * i + j] = min((i + (rng.RandomFloat())) * delta,
                                        OneMinusEpsilon);

    // Permute LHS samples in each dimension
    for (uint32_t i = 0; i < nDim; ++i) {
        for (uint32_t j = 0; j < nSamples; ++j) {
            uint32_t other = j + (rng.RandomUInt() % (nSamples - j));
            swap(samples[nDim * j + i], samples[nDim * other + i]);
        }
    }
}

// multiple cores and threading
static pthread_t *threads;
static Mutex *taskQueueMutex = Mutex::Create();
static std::vector<Task *> taskQueue;
static Semaphore *workerSemaphore;
static uint32_t numUnfinishedTasks;
static ConditionVariable *tasksRunningCondition;
static void *taskEntryPoint(void *arg);

// Mutex
Mutex *Mutex::Create()
{
    int sz = sizeof(Mutex);
    sz = (sz + (L1_CACHE_LINE_SIZE-1)) & ~(L1_CACHE_LINE_SIZE-1);
    return new (AllocAligned(sz)) Mutex;
}

void Mutex::Destroy(Mutex *m)
{
    m->~Mutex();
    FreeAligned(m);
}

Mutex::Mutex()
{
    if (pthread_mutex_init(&mutex, NULL)) {
        cout << "ERROR: pthread_mutex_init." << endl;
				exit(1144);
		}
}

Mutex::~Mutex()
{
    if (pthread_mutex_destroy(&mutex)) {
        cout << "ERROR: pthread_mutex_destroy." << endl;
				exit(1142);
		}
}

// MutexLock
MutexLock::MutexLock(Mutex &m) : mutex(m)
{
    if (pthread_mutex_lock(&m.mutex)) {
        cout << "ERROR: pthread_mutex_lock." << endl;
				exit(1143);
		}
}

MutexLock::~MutexLock()
{
    if (pthread_mutex_unlock(&mutex.mutex)) {
        cout << "ERROR: pthread_mutex_unlock." << endl;
				exit(1145);
		}
}

// Semaphore
Semaphore::Semaphore()
{
    char name[32];
    sprintf(name, "areport.%d-%d", (int)getpid(), count++);
    sem = sem_open(name, O_CREAT, S_IRUSR|S_IWUSR, 0);
    if (!sem) {
        cout << "ERROR: sem_open " << sem << endl;
				exit(1146);
		}
}

int Semaphore::count = 0;

Semaphore::~Semaphore()
{
    if (sem_close(sem)) {
        cout << "ERROR: sem_close." << endl;
				exit(1146);
		}
}

void Semaphore::Wait()
{
    if (sem_wait(sem)) {
        cout << "ERROR: sem_wait." << endl;
				exit(1147);
		}
}

bool Semaphore::TryWait()
{
    return (sem_trywait(sem) == 0);
}

void Semaphore::Post(int count)
{
    while (count-- > 0)
		{
        if (sem_post(sem)) {
						cout << "ERROR: sem_post." << endl;
						exit(1148);
				}
		}
}

// ConditionVariable
ConditionVariable::ConditionVariable()
{
   int err;
   if ((err = pthread_cond_init(&cond, NULL)) != 0) {
        cout << "ERROR: cv pthread_cond_init " << err << endl;
				exit(1154);
		}
   if ((err = pthread_mutex_init(&mutex, NULL)) != 0) {
        cout << "ERROR: cv pthread_mutex_init " << err << endl;
				exit(1153);
		}
}

ConditionVariable::~ConditionVariable()
{
    pthread_cond_destroy(&cond);
    pthread_mutex_destroy(&mutex);
}

void ConditionVariable::Lock()
{
    int err;
    if ((err = pthread_mutex_lock(&mutex)) != 0) {
        cout << "ERROR: cv_lock " << err << endl;
				exit(1149);
		}
}

void ConditionVariable::Wait()
{
    int err;
    if ((err = pthread_cond_wait(&cond, &mutex)) != 0) {
        cout << "ERROR: cv_wait " << err << endl;
				exit(1150);
		}
}

void ConditionVariable::Signal()
{
    int err;
    if ((err = pthread_cond_signal(&cond)) != 0) {
        cout << "ERROR: cv_signal " << err << endl;
				exit(1151);
		}
}

void ConditionVariable::Unlock() {
    int err;
    if ((err = pthread_mutex_unlock(&mutex)) != 0) {
        cout << "ERROR: cv_unlock " << err << endl;
				exit(1152);
		}
}

// thread entry point
static void *taskEntryPoint(void *arg)
{
		IF_DEBUG(cout << "taskEntyPoint()" << endl);
		
    while (true)
		{
        workerSemaphore->Wait();
				
        // get task from task queue
        Task *myTask = NULL;
        { MutexLock lock(*taskQueueMutex);
				
        if (taskQueue.size() == 0)
            break;
						
        myTask = taskQueue.back();
        taskQueue.pop_back();
        }
				
				IF_DEBUG(cout << " Acquired Task, Running." << endl);

        // run acquired rendererTask
        myTask->Run();
        tasksRunningCondition->Lock();
				
        int unfinished = --numUnfinishedTasks;
        if (unfinished == 0)
            tasksRunningCondition->Signal();
						
        tasksRunningCondition->Unlock();
    }
		
    // Cleanup from task thread and exit
    pthread_exit(NULL);
    return 0;
}

void startTasks(const vector<Task *> &tasks)
{
		// sequential execution requested
		if (Config.nCores == 1) {
				for (int i=0; i < tasks.size(); i++) {
						IF_DEBUG(cout << "startTasks:: Starting SEQUENTIAL run task [" << i << "]" << endl);
						tasks[i]->Run();
				}
				return;
		}
		
		// init
		if (!threads)
				TasksInit();
				
		// update taskQueue
		{ MutexLock lock(*taskQueueMutex);
		for (unsigned int i = 0; i < tasks.size(); ++i)
				taskQueue.push_back(tasks[i]);
		}
		
		// update running flag
		tasksRunningCondition->Lock();
		numUnfinishedTasks += tasks.size();
		tasksRunningCondition->Unlock();

		// update counter
		workerSemaphore->Post(tasks.size());
}

void waitUntilAllTasksDone()
{
		// nothing to do if sequential execution requested
		if (Config.nCores == 1)
				return;
				
		// make sure TasksInit() was called
		if (!tasksRunningCondition)
				return;
				
		// update running flag
		tasksRunningCondition->Lock();
		while (numUnfinishedTasks > 0)
				tasksRunningCondition->Wait();
		tasksRunningCondition->Unlock();
}

Task::~Task() {
}

void TasksInit()
{
		if (Config.nCores == 1)
				return;

		// setup semaphore and check
    static const int nThreads = numberOfCores();
    workerSemaphore = new Semaphore;
    tasksRunningCondition = new ConditionVariable;				
				
		// create threads
    threads = new pthread_t[nThreads];
    for (int i = 0; i < nThreads; ++i) {
				IF_DEBUG(cout << "TasksInit:: Creating new thread [" << i << "]" << endl);
        int err = pthread_create(&threads[i], NULL, &taskEntryPoint, reinterpret_cast<void *>(i));
        if (err != 0) {
            cout << "Renderer::TasksInit(): ERROR from pthread_create [" << err << "]" << endl;
						exit(1141);
				}
    }
}

void TasksCleanup()
{
		if (Config.nCores == 1)
				return;
				
    if (!taskQueueMutex || !workerSemaphore)
        return;
				
		// verify task queue is empty
    { MutexLock lock(*taskQueueMutex);
    if(taskQueue.size() != 0)
				exit(1155);
    }

		// increment worker counter
    static const int nThreads = numberOfCores();
		
    if (workerSemaphore != NULL)
        workerSemaphore->Post(nThreads);

		// close handles and delete threads
    if (threads != NULL) {
        for (int i = 0; i < nThreads; ++i) {
						IF_DEBUG(cout << "TasksCleanup:: Joining thread [" << i << "]" << endl);
            int err = pthread_join(threads[i], NULL);
            if (err != 0) {
                cout << "ERROR: tasks_cleanup pthread_join " << err << endl;
								exit(1156);
						}
        }
        delete[] threads;
        threads = NULL;
    }
}

int numberOfCores()
{
    if (Config.nCores > 0)
				return Config.nCores;

    return sysconf(_SC_NPROCESSORS_ONLN);
}
