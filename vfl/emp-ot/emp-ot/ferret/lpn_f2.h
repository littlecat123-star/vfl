#ifndef EMP_LPN_F2K_H__
#define EMP_LPN_F2K_H__

#include "emp-tool/emp-tool.h"
using namespace emp;

//Implementation of local linear code on F_2^k
//Performance highly dependent on the CPU cache size
template<typename IO, int d = 10>
class LpnF2 { public:
	int party;
	int k, n;
	ThreadPool * pool;
	IO *io;
	int threads;
	block seed;
	int mask;
	LpnF2 (int party, int n, int k, ThreadPool * pool, IO *io, int threads) {
		this->party = party;
		this->k = k;
		this->n = n;
		this->pool = pool;
		this->io = io;
		this->threads = threads;
		mask = 1;
		while(mask < k) {
			mask <<=1;
			mask = mask | 0x1;
		}
	}

	void __compute4(block * nn, const block * kk, int i, PRP * prp) {
		block tmp[10];
		for(int m = 0; m < 10; ++m)
			tmp[m] = makeBlock(i, m);
		AES_ecb_encrypt_blks(tmp, 10, &prp->aes);
		uint32_t* r = (uint32_t*)(tmp);
		for(int m = 0; m < 4; ++m)
			for (int j = 0; j < d; ++j) {
				int index = (*r) & mask;
				++r;
				index = index >= k? index-k:index;
				nn[i+m] = nn[i+m] ^ kk[index];
			}
	}

	void __compute1(block * nn, const block * kk, int i, PRP*prp) {
		block tmp[3];
		for(int m = 0; m < 3; ++m)
			tmp[m] = makeBlock(i, m);
		prp->permute_block(tmp, 3);
		uint32_t* r = (uint32_t*)(tmp);
		for (int j = 0; j < d; ++j)
			nn[i] = nn[i] ^ kk[r[j]%k];
	}

	void task(block * nn, const block * kk, int start, int end) {
		PRP prp(seed);
		int j = start;
		for(; j < end-4; j+=4)
			__compute4(nn, kk, j, &prp);
		for(; j < end; ++j)
			__compute1(nn, kk, j, &prp);
	}

	void compute(block * nn, const block * kk) {
		vector<std::future<void>> fut;
		int width = n/threads;
		seed = seed_gen();
		for(int i = 0; i < threads - 1; ++i) {
			int start = i * width;
			int end = min((i+1)* width, n);
			fut.push_back(pool->enqueue([this, nn, kk, start, end]() {
				task(nn, kk, start, end);
			}));
		}
		int start = (threads - 1) * width;
		int end = min(threads * width, n);
		task(nn, kk, start, end);

		for (auto &f: fut) f.get();
	}

	block seed_gen() {
		block seed;
		if(party == ALICE) {
			PRG prg;
			prg.random_block(&seed, 1);
			io->send_data(&seed, sizeof(block));
		} else {
			io->recv_data(&seed, sizeof(block));
		}io->flush();
		return seed;
	}
	void bench(block * nn, const block * kk) {
		vector<std::future<void>> fut;
		int width = n/threads;
		for(int i = 0; i < threads - 1; ++i) {
			int start = i * width;
			int end = min((i+1)* width, n);
			fut.push_back(pool->enqueue([this, nn, kk, start, end]() {
				task(nn, kk, start, end);
			}));
		}
		int start = (threads - 1) * width;
		int end = min(threads * width, n);
		task(nn, kk, start, end);

		for (auto &f: fut) f.get();
	}

};
#endif
