//This realizes the paper's online phase
#include "ttp.h"
#include "hash.h"
#include "inv.h"
#include "PackedArray.h"
#include <omp.h>

#define kk 4  // new,as size distance between A and C

using namespace emp;
using namespace std;
int port, party;

int main(int argc, char** argv) {
  omp_set_num_threads(THREADS);
  
    if (argc!=3) {
        cout <<"Must supply party and port number"<<endl;
    }

    parse_party_and_port(argv, &party, &port);

    cout <<"Buckets: "<<BUCKETS<<", BETA: "<<BETA<<", #Hash elements: "<<HASHELEMENTS<<", modulus: "<<MODULUS<<endl;
    
    
    //Only known to Alice
    BT *sA = (BT *) calloc(1, BUCKETS*sizeof(BT));
    BT *rA = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT));
    PackedArray* toBobPacked = NULL;  
    PackedArray* toAndyPacked = NULL; // NEW,send the finnal payload
    BT *aliceCuckooTable = NULL;
    BT *toBob = NULL;
    BT *fromBob = NULL;
    PackedArray* fromBobPacked = NULL;
    //new
    BT *fromBobPayload = NULL;
    PackedArray* fromBobPayloadPacked = NULL;
    BT *aliceTable2 = NULL; //save the finall payload table
    // before two item are used to save the payload from bob

    //let BOb be the corrordiantor (NEW)
    BT *bobTable2 = NULL;//save the 2-th ID hash
    BT *AndyTable = NULL; //save the andy's final payload table
    // end new

    //Only known to Bob 
    BT *INV = (BT*) calloc(1, MODULUS*sizeof(BT));
    BT *bobTable = NULL; // save the 1-th ID hash
    BT *rB = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT));
    BT *sB = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT));
    PackedArray* fromAlicePacked = NULL;
    
    BT *fromAlice = NULL;
    BT *fromAlicePayload = NULL;
    PackedArray* toAlicePacked = NULL;
    BT *toAlicePayload = NULL;
    PackedArray* toAlicePayloadPacked = NULL; //save corresponding payload

    auto quad = clock_start();
    if (party == BOB) {
      computeInverses(INV);
      cout <<"Computing inverses in "<<time_from(quad)/1000<<" ms"<<endl;
    }
    //For simplicity, we quickly generate the same triples on the fly on both parties. This is typically done offline before.
    cout <<"Computing triples "<<flush;
    quad = clock_start();
    computeTriples(rA, sA, rB, sB);
    cout <<"in "<<time_from(quad)/1000<<" ms"<<endl;

    if (party == ALICE) {
      free(rB);
      free(sB);
      free(INV);
    } else {
      free(rA);
      free(sA);
      }
       
    //Fix seed for debugging of fixed input
    block seed = makeBlock(7, 8);
    PRG prg;
    prg.reseed(&seed);

    //We use Cuckoo hashing with 3 hashfunctions. We can prepare these offline.
    BT *h1 = (BT *) calloc(1, HASHELEMENTS*sizeof(BT));
    BT *h2 = (BT *) calloc(1, HASHELEMENTS*sizeof(BT));
    BT *h3 = (BT *) calloc(1, HASHELEMENTS*sizeof(BT));
    
    cout <<"Preparing hash functions..."<<endl;
    auto hashFTime = clock_start();
    computeHashFunctions(prg, h1, h2, h3);
    cout <<"Time to prepare hash functions: "<<time_from(hashFTime)/1000<<" ms"<<endl;
    
    double hashing = 0.0;
    double CPUTime = 0.0;
    
    if(party == ALICE) {
      //Now we start preparing the Cuckoo table
        toBobPacked = PackedArray_create(BITS, BUCKETS);
        aliceCuckooTable = (BT *) calloc(1, BUCKETS*sizeof(BT));
        toBob = (BT *) calloc(1, BUCKETS*sizeof(BT));
        BT *aliceInput  = (BT*) calloc(1, (NELEMENTS/kk) * sizeof(BT));
        fromBob = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT));
        fromBobPacked = PackedArray_create(BITS, BETA*BUCKETS);
        fromBobPayload = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT)); // NEW
        fromBobPayloadPacked = PackedArray_create(BITS, BETA*BUCKETS); // NEW
        aliceTable2 = (BT *) calloc(1, 3*BUCKETS*sizeof(BT)); // new finally table size,3 hash functions
        //Fill Alice with fixed input
	for (size_t i = 0; i<NELEMENTS/kk; i++) {
	  aliceInput[i] = i;
	}
        cout<<"NELEMENTS="<<NELEMENTS/kk<<endl;
	
        cuckooEntry **bufAliceCuckooTable = (cuckooEntry **) calloc(1, BUCKETS*sizeof(cuckooEntry*));

        auto start = clock_start();

	//Populate Cuckoo table...
        computeCuckoo(aliceInput, h1, h2, h3, bufAliceCuckooTable);

        //...and convert it into a proper format to send it
        for (size_t i = 0; i < BUCKETS; i++) {
            if (bufAliceCuckooTable[i]!=NULL) {
                aliceCuckooTable[i] = (bufAliceCuckooTable[i])->value;
            } else {
	      //Mark as invalid
                aliceCuckooTable[i] = 3*HASHELEMENTS;
            }
        } // fill the cuckoo hash table

	hashing = time_from(start);
	CPUTime += hashing;
        cout <<"Alice Cuckoo hashing time: "<<hashing/1000<<" ms"<<endl;

	//No use for hash functions anymore
	free(h1);
	free(h2);
	free(h3);

	free(bufAliceCuckooTable);
	
        //Compute $c$s for Bob
	double elapsed = 0.0;
        // double compute_time = 0.0;
        start=clock_start();
        for (size_t i = 0; i < BUCKETS; i++) {
            toBob[i] = (sA[i] - aliceCuckooTable[i] + MODULUS) % MODULUS;
        }
        PackedArray_pack(toBobPacked, 0, toBob, BUCKETS);
	elapsed = time_from(start);
	CPUTime += elapsed;
        cout <<"Alice compute and pack c time: "<<elapsed/(1000)<<" ms"<<endl;

	free(aliceCuckooTable);
	free(aliceInput);
	
        //END ALICE
    } else {
        //BOB

      //We prepare a regular hash table
        bobTable = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT));
        // bobTable2 = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT)); // init the second hash table
        fromAlicePacked = PackedArray_create(BITS, BUCKETS);
        // PackedArray* fromAlicePayloadPacked = NULL; // new as Andy
        // fromAlicePayloadPacked = PackedArray_create(BITS, 3*BUCKETS);//new
        fromAlice = (BT *) calloc(1, BUCKETS*sizeof(BT));
        fromAlicePayload = (BT *) calloc(1, 3*BUCKETS*sizeof(BT));
        toAlicePacked = PackedArray_create(BITS, BETA*BUCKETS);
        toAlicePayload = (BT *) calloc(1, BETA*BUCKETS*sizeof(BT)); // new
	toAlicePayloadPacked = PackedArray_create(BITS, BETA*BUCKETS); // save the payload
        //Compute fixed input
        BT *bobInput  = (BT*) malloc(NELEMENTS * sizeof(BT));
        
        for (size_t i = 0; i<NELEMENTS; i++) {
	  bobInput[i] = i;
	}  // fill with Bob's input
	cout<<"NELEMENTS="<<NELEMENTS<<endl;
        bucket *_bobTable = (bucket*) calloc(1, BUCKETS*sizeof(bucket));
        bucket *_toAlicePayload = (bucket*) calloc(1, BUCKETS*sizeof(bucket)); // new
        // bucket *_bobTable2 = (bucket*) calloc(1, BUCKETS*sizeof(bucket)); // new
        auto t = clock_start();

	//Compute regular hash table
	computeHash(bobInput, h1, h2, h3, _bobTable, bobTable,_toAlicePayload, toAlicePayload); // bobTable is the simple hash table
	// computeHash(bobInput, h1, h2, h3, _bobTable2, bobTable2,); // new
	hashing = time_from(t);
	CPUTime += hashing;
        cout <<"Bob hashing and computing payload: "<<hashing/1000<<" ms"<<endl; // double time
    
	//No use for hash functions anymore
	free(h1);
	free(h2);
	free(h3);

	free(_bobTable);
        free(_toAlicePayload);
        // free(_bobTable2);
	free(bobInput);
        //END BOB
    }

    cout <<"Connecting..."<<endl;
    double total_time = 0.0;  // only consider the time when alice and bob begin to interact
    double compute_time = 0.0;

    auto ttime_including_idle = clock_start();  // waiting time

    //Server or local environment: NetIO* io = new NetIO(party==ALICE?nullptr:"127.0.0.1", port);
    //West Coast: NetIO* io = new NetIO(party==ALICE?nullptr:"104.42.191.30", port);
    //EU: NetIO* io = new NetIO(party==ALICE?nullptr:"13.95.162.209", port);

    
    NetIO* io = new NetIO(party==ALICE?nullptr:"127.0.0.1", port);
    
    if (party == ALICE) {

      //Send the $c$s
      
        //4 = number of bytes per uint32_t slot
        size_t myLength = 4*PackedArray_bufferSize(toBobPacked);
        size_t myLength2 = 4*PackedArray_bufferSize(toBobPacked); // new
        cout <<"Alice sending "<<myLength<<" Bytes"<<endl;
        io->sync();
	
	ttime_including_idle = clock_start();
        auto start=clock_start();
        mySend(io, (char*)toBobPacked, myLength);
        auto elapsed = time_from(start);
        total_time +=elapsed;
        cout <<"Alice sent c in "<<elapsed/1000<<" ms"<<endl;

	free(toBobPacked);
        
	//Alice gets back $d$s
        myLength = 4*PackedArray_bufferSize(fromBobPacked);
        myLength2 = 4*PackedArray_bufferSize(fromBobPayloadPacked);
        
        cout <<"Alice receiving "<<myLength+myLength2<<" Bytes"<<endl; // double length
        io->sync();     
        start=clock_start();
        myRecv(io, (char*) fromBobPacked, myLength);       
	auto pack = clock_start();
        PackedArray_unpack(fromBobPacked, 0, fromBob, BETA * BUCKETS);
        CPUTime += time_from(pack);	
        elapsed = time_from(start);
        total_time += elapsed;
        // cout<<"Alice receive d"<<endl;
        cout <<"Alice receive and unpack d  in "<<elapsed/1000<<" ms"<<endl;
        io->sync();     
        start=clock_start();
        myRecv(io, (char*) fromBobPayloadPacked, myLength2); //new
        auto packpayload = clock_start();
        PackedArray_unpack(fromBobPayloadPacked, 0, fromBobPayload, BETA * BUCKETS);
	CPUTime += time_from(packpayload);	
        elapsed = time_from(start);
        total_time +=elapsed;
        // cout<<"Alice receive payload"<<endl;
        cout <<"Alice receive and unpack payload in "<<elapsed/1000<<" ms"<<endl;
	free(fromBobPacked);
        free(fromBobPayloadPacked); // new
        
        size_t match = 0;
        cout<<"Alice begin to compare"<<endl;
        //Compute intersection
        start=clock_start();
	#pragma omp parallel for reduction(+:match)
        for (size_t i = 0; i < BUCKETS; i++) {
	  for (size_t j = 0; j < BETA; j++) {
                if (rA[i*BETA + j]==fromBob[i * BETA + j]) {
                    match++;
                    // new
                    auto startcompu = clock_start();
                    for (size_t k = 0; k < 3; k++) 
                    {
		       aliceTable2[i*3 + k]=fromBobPayload[i*BETA + j]+k ; 
                    } 
                    auto end_compute=time_from(startcompu);
                    compute_time +=end_compute;  // need to record compute time          
                    }
            }
        }
        elapsed = time_from(start);
        total_time +=elapsed;
	CPUTime+=elapsed;
        cout <<"Alice compare r_A and d in "<<(elapsed-compute_time)/1000<<" ms, size "<<match<<endl;
        cout <<"A<-->C:CPUtime="<<(CPUTime-compute_time)/1000<<" ms,total_time="<<(total_time-compute_time)/1000<<" ms"<<endl;
        cout <<"A<-->C:total time(including waiting+hashing)="<<(time_from(ttime_including_idle)-compute_time+hashing)/1000<<" ms"<<endl;
        cout <<"A<-->C finish, A<---->B begin"<<endl;
        cout <<"A<-->B, A compute final payload time: "<<compute_time/1000<<" ms"<<endl;
        
        toAndyPacked = PackedArray_create(BITS, BUCKETS*3); // NELEMENTS/kk  
        pack = clock_start();
        PackedArray_pack(toAndyPacked, 0, aliceTable2, BUCKETS);
        elapsed = time_from(pack);
        compute_time +=elapsed;
        total_time +=elapsed;
        CPUTime +=elapsed;
        cout <<"A<-->B CPUTime="<< compute_time/1000<<" ms"<<endl;
	// send the finnal payload to compare  
        myLength = 4*PackedArray_bufferSize(toAndyPacked); 
        cout <<"Alice sending the finnal payload to party B "<<myLength<<" Bytes"<<endl;      
        io->sync();
        start=clock_start();
	mySend(io, (char*) toAndyPacked, myLength);	
        elapsed = time_from(start);
        compute_time += elapsed;
        total_time +=elapsed;
        // compute_time +=elapsed;
        cout <<"Alice send finall payload in "<<elapsed/1000<<" ms"<<endl;
        cout <<"A<-->B total time="<< compute_time/1000<<" ms"<<endl;
	free(toAndyPacked);
        
	//END ALICE

    } else {
        //BOB

      //Bob gets the $c$s
      
        //4 = number of bytes per uint32_t slot
        size_t myLength = 4*PackedArray_bufferSize(fromAlicePacked);
        size_t myLength2 = myLength; // new 
        //BOB
        cout <<"Bob receiving "<<myLength<<" Bytes"<<endl;
        io->sync();
	
	ttime_including_idle = clock_start();
        auto start=clock_start();

        //io->recv_data(fromAlicePacked, myLength);
	myRecv(io, (char *) fromAlicePacked, myLength); // cannot receive?
        //cout<<"receiving c"<<endl;
	auto unpack = clock_start();
        PackedArray_unpack(fromAlicePacked, 0, fromAlice, BUCKETS);
        // cout<<"unpack the c"<<endl;
	auto elapsed = time_from(start);
	CPUTime +=time_from(unpack);
        total_time +=elapsed;
        cout <<"Bob receive and unpack c in "<<elapsed/1000<<" ms"<<endl;

	free(fromAlicePacked);
        
        //Compute result: d =  (fromAlice + myHash + sB) / rB 
        start = clock_start();
	bobsOperations(bobTable, fromAlice, rB, sB, INV); // compute d,d is saved in bobTable
	PackedArray_pack(toAlicePacked, 0, bobTable, BETA*BUCKETS); // pack
        // cout<<"BOb pack d"<<endl;
	PackedArray_pack(toAlicePayloadPacked, 0, toAlicePayload, BETA*BUCKETS); // new
        // cout<<"BOb pack payload"<<endl;
        elapsed = time_from(start);
        total_time +=elapsed;
	CPUTime+=elapsed;
        cout <<"Bob computes and pack d,payload in "<<elapsed/1000<<" ms"<<endl;
        
	free(fromAlice);
	free(bobTable);      
        free(toAlicePayload); // NEW


	//Send all $d$s back
        myLength = 4*PackedArray_bufferSize(toAlicePacked);
        myLength2 = 4*PackedArray_bufferSize(toAlicePayloadPacked); // new
        cout <<"Bob sending "<<myLength+myLength2<<" Bytes"<<endl;
        io->sync();
        start=clock_start();
	mySend(io, (char*) toAlicePacked, myLength);	
        elapsed = time_from(start);
        total_time +=elapsed;
        cout <<"Bob sent d in "<<elapsed/1000<<" ms"<<endl;
	free(toAlicePacked);
        io->sync();
        start=clock_start();
	mySend(io, (char*) toAlicePayloadPacked, myLength2);	
        elapsed = time_from(start);
        total_time +=elapsed;
        cout <<"Bob sent payload in "<<elapsed/1000<<" ms"<<endl;       
        free(toAlicePayloadPacked);
        cout <<"A<-->C:CPUtime="<<CPUTime/1000<<" ms,total_time="<<total_time/1000<<" ms"<<endl;
        cout <<"A<-->C:total time(including waiting+hashing)="<<(time_from(ttime_including_idle)+hashing)/1000<<" ms"<<endl;
        
        cout <<"below,Bob serve as Andy"<<endl;
        // new 
        AndyTable = (BT *) calloc(1, 3*BUCKETS*sizeof(BT)); 
        PackedArray* fromAlicePayloadPacked = NULL; // new as Andy
        fromAlicePayloadPacked = PackedArray_create(BITS, 3*BUCKETS);//new    
        myLength = 4*PackedArray_bufferSize(fromAlicePayloadPacked);
        // cout <<"Andy receiving "<<myLength<<" Bytes"<<endl;
        start=clock_start();
        io->sync();  
        //io->recv_data(fromAlicePacked, myLength);
	myRecv(io, (char *) fromAlicePayloadPacked, myLength); 
        // cout<<"Andy receiving payload from alice "<<endl;
        unpack = clock_start();
        PackedArray_unpack(fromAlicePayloadPacked, 0, fromAlicePayload, 3*BUCKETS);
        CPUTime +=time_from(unpack);
        cout <<"Andy unpack payload in "<<time_from(unpack)/1000<<" ms"<<endl;      
        elapsed = time_from(start);   
        compute_time +=elapsed;
        free(fromAlicePayloadPacked); 
        cout<<"Andy receive final and unpack payload  in "<<elapsed/1000<<" ms"<<endl;
 
        for (size_t i = 0; i < BUCKETS; i++) {
	  for (size_t j = 0; j < 3; j++) {
                AndyTable[i*3+j] = fromAlicePayload[i*3+j]+9;
          }
        } // init andy's local final payload
        size_t Andy_match = 0;
        size_t flag = 0;
        start = clock_start();
        #pragma omp parallel for reduction(+:Andy_match)
        for (size_t i = 0; i < BUCKETS; i++) {
	  for (size_t j = 0; j < 3; j++) {
                for (size_t ll = 0; ll < 3*BUCKETS; ll++) {
                if (AndyTable[i*3+ j]==fromAlicePayload[ll]) {
                    Andy_match++;
                    flag = 1;
                }
                if (flag==1){ 
                break;}
        }          
          if (flag==1){ 
                flag=0;
                break;}
        }
        }
        elapsed = time_from(start);
        compute_time +=elapsed;
        CPUTime +=elapsed;
        cout<<"Andy compare locally in "<<elapsed/1000<<" ms, match:"<<Andy_match<<endl;
        cout<<"A<-->B:Andy total time:"<<compute_time/1000<<" ms"<<endl;
    }

    cout <<"Total time for party "<<myName<<": "<<total_time/1000<<" ms, sent: "<<(io->counter)<<" Byte, i.e., "<<((double)io->counter)/(1024*1024)<<" MByte"<<endl;
    // running time after hashing without waiting
    cout <<"Total time including waiting for "<<myName<<": "<<time_from(ttime_including_idle)/1000<<" ms"<<endl;
    // running time after hashing with waiting
    cout <<"CPU time for "<<myName<<": "<<(CPUTime)/1000<<" ms"<<endl;
    // computing time 
    // cout<<"hashing time "<<hashing/1000<<"ms"<<endl;
    cout <<"Total time including waiting and hashing for "<<myName<<": "<<(time_from(ttime_including_idle)+hashing)/1000<<" ms"<<endl;
    // hash time+running time
    cout <<"---------"<<endl;
}