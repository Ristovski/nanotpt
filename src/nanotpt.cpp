/**
	This file is part of The Powder Toy.
	The Powder Toy is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	The Powder Toy is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	You should have received a copy of the GNU General Public License
	along with The Powder Toy.  If not, see <https://www.gnu.org/licenses/>.
**/

#include <algorithm>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

#define SIMULATIONW 800
#define SIMULATIONH 600

#define TYPE_NONE 0
#define TYPE_SOLID 1
#define TYPE_POWDER 2
#define TYPE_LIQUID 3
#define TYPE_GAS 4
#define TYPE_PARTICLE 5

#define PART(x, y) (x) + ((y) * SIMULATIONW)

#define PART_POS_QUANT(x) ((int)(x + 0.5f))

bool displacementMatrix[6][6]{
	{false, false, false, false, false, false},
	{false, false, false, false, false, false},
	{true, false, false, true, true, false},
	{true, false, false, false, true, false},
	{true, false, false, false, true, false},
	{true, true, true, true, true, true}};

struct atom
{
	uint8_t type = TYPE_NONE;
	float vx = 0.0f;
	float vy = 0.0f;
	float x = 0.0f;
	float y = 0.0f;
	bool mutex = false;
};

struct region_bounds
{
	int x;
	int y;
	int w;
	int h;
};

#define GRAVITYAY 0.5f
#define VLOSS 0.99f
#define DIFFUSION 0.2f

#define ISTP 1
#define COLLISIONLOSS 0.1f

float randfd()
{
	return ((rand() % 1000) / 500.0f) - 1.0f;
}

int randd()
{
	return (rand() % 2) * 2 - 1;
}

bool do_move(atom *parts, atom &current, float resultx, float resulty)
{
	int resultx_quant = PART_POS_QUANT(resultx);
	int resulty_quant = PART_POS_QUANT(resulty);
	if (resultx_quant < 0 || resultx_quant >= SIMULATIONW || resulty_quant < 0 || resulty_quant >= SIMULATIONH)
	{
		current.type = TYPE_NONE;
		return true;
	}

	atom &target = parts[PART(resultx_quant, resulty_quant)];

	if (displacementMatrix[current.type][target.type])
	{
		atom temp = target;
		target = current;
		current = temp;
		target.x = resultx;
		target.y = resulty;
		return true;
	}
	else
	{
		return false;
	}
}

std::atomic<uint32_t> last_partcount(0);

void simulate_region(atom *parts, region_bounds region, bool mutex)
{
	int nx, ny, neighbourSpace, neighbourDiverse;
	bool neighbourBlocking;

	float mv = 0.0f, resultx = 0.0f, resulty = 0.0f;
	int resultx_quant, resulty_quant;

	for (int gridY = region.y; gridY < region.y + region.h; gridY++)
	{
		if (gridY == 0 || gridY == SIMULATIONH - 1)
			continue;
		for (int gridX = region.x; gridX < region.x + region.w; gridX++)
		{
			if (gridX == 0 || gridX == SIMULATIONW - 1)
				continue;

			atom &current = parts[PART(gridX, gridY)];

			if (current.type == TYPE_NONE)
				continue;

			last_partcount++;

			if (current.mutex == mutex)
				continue;

			current.mutex = mutex;

			if (current.type == TYPE_GAS || current.type == TYPE_POWDER || current.type == TYPE_LIQUID)
			{
				current.vx *= VLOSS;
				current.vy *= VLOSS;
			}

			if (current.type == TYPE_POWDER || current.type == TYPE_LIQUID)
			{
				current.vy += GRAVITYAY;
			}

			if (current.type == TYPE_GAS)
			{
				current.vx += DIFFUSION * randfd();
				current.vy += DIFFUSION * randfd();
			}

			if (current.type == TYPE_LIQUID)
			{
				current.vx += DIFFUSION * randfd() * 0.1f;
				current.vy += DIFFUSION * randfd() * 0.1f;
			}

			neighbourSpace = neighbourDiverse = 0;
			neighbourBlocking = true;

			for (nx = -1; nx < 2; nx++)
			{
				for (ny = -1; ny < 2; ny++)
				{
					if (nx || ny)
					{
						atom &neighbour = parts[PART(gridX + nx, gridY + ny)];
						if (neighbour.type == TYPE_NONE)
						{
							neighbourSpace++;
							neighbourBlocking = false;
						}
						if (neighbour.type != current.type)
							neighbourDiverse++;
						if (displacementMatrix[neighbour.type][current.type])
							neighbourBlocking = false;
					}
				}
			}

			if (neighbourBlocking)
			{
				current.vx = 0.0f;
				current.vy = 0.0f;
				continue;
			}

			if ((fabsf(current.vx) <= 0.01f && fabsf(current.vy) <= 0.01f) || current.type == TYPE_SOLID)
				continue;

			mv = fmaxf(fabsf(current.vx), fabsf(current.vy));

			//if (mv < ISTP)
			{
				resultx = current.x + current.vx;
				resulty = current.y + current.vy;
			}
			//else
			{
				//Interpolation, TODO
			}

			resultx_quant = PART_POS_QUANT(resultx);
			resulty_quant = PART_POS_QUANT(resulty);

			int clearx = gridX;
			int cleary = gridY;

			float clearxf = current.x;
			float clearyf = current.y;

			if (resultx_quant != gridX || resulty_quant != gridY)
			{
				if (do_move(parts, current, resultx, resulty))
					continue;
				if (current.type == TYPE_GAS)
				{
					if (do_move(parts, current, 0.25f + (float)(2 * gridX - resultx_quant), 0.25f + resulty_quant))
					{
						current.vx *= COLLISIONLOSS;
						continue;
					}
					else if (do_move(parts, current, 0.25f + resultx_quant, 0.25f + (float)(2 * gridY - resulty_quant)))
					{
						current.vy *= COLLISIONLOSS;
						continue;
					}
					else
					{
						current.vx *= COLLISIONLOSS;
						current.vy *= COLLISIONLOSS;
						continue;
					}
				}
				if (current.type == TYPE_LIQUID || current.type == TYPE_POWDER)
				{
					if (resultx_quant != gridX && do_move(parts, current, resultx, gridY))
					{
						current.vx *= COLLISIONLOSS;
						current.vy *= COLLISIONLOSS;
						continue;
					}
					else if (resulty_quant != gridY && do_move(parts, current, gridX, resulty))
					{
						current.vx *= COLLISIONLOSS;
						current.vy *= COLLISIONLOSS;
						continue;
					}
					else
					{
						int scanDirection = randd();
						if (clearx != gridX || cleary != gridY || neighbourDiverse || neighbourSpace)
						{
							float dx = current.vx - current.vy * scanDirection;
							float dy = current.vy + current.vx * scanDirection;
							if (fabsf(dy) > fabsf(dx))
								mv = fabsf(dy);
							else
								mv = fabsf(dx);
							dx /= mv;
							dy /= mv;
							if (do_move(parts, current, clearxf + dx, clearyf + dy))
							{
								current.vx *= COLLISIONLOSS;
								current.vy *= COLLISIONLOSS;
								continue;
							}
							float swappage = dx;
							dx = dy * scanDirection;
							dy = -swappage * scanDirection;
							if (do_move(parts, current, clearxf + dx, clearyf + dy))
							{
								current.vx *= COLLISIONLOSS;
								current.vy *= COLLISIONLOSS;
								continue;
							}
						}
						current.vx *= COLLISIONLOSS;
						current.vy *= COLLISIONLOSS;
					}
				}
			}
		}
	}
}

int threadcount = 0;
int regioncount = 0;
int region_group_count = 0;
region_bounds *regions;
region_bounds **region_groups;
region_bounds *active_regions;
std::thread *threads;
std::atomic_flag **locks;

bool mutex = true;

bool exiting = false;

std::atomic<uint8_t> runs;

void simulate_region_thread(std::atomic_flag *lock, atom *parts, int threadid)
{
	//while (lock->test_and_set(std::memory_order_acquire));
	while (true)
	{
		while (lock->test_and_set(std::memory_order_acquire))
			;
		if (exiting)
			break;
		simulate_region(parts, active_regions[threadid], mutex);
		runs++;
	}
}

void init_simulation(int threadcount_, int groupcount_, atom *parts)
{
	threadcount = threadcount_;
	region_group_count = std::min(groupcount_, threadcount);
	regioncount = threadcount_ * region_group_count;

	locks = new std::atomic_flag *[threadcount];
	threads = new std::thread[threadcount];
	for (int i = 0; i < threadcount; i++)
	{
		locks[i] = new std::atomic_flag();
		locks[i]->test_and_set(std::memory_order_acquire);
		threads[i] = std::thread(simulate_region_thread, locks[i], parts, i);
	}

	regions = new region_bounds[regioncount];
	region_groups = new region_bounds *[region_group_count];
	for (int i = 0; i < region_group_count; i++)
		region_groups[i] = new region_bounds[threadcount];

	int regionwidth = SIMULATIONW / regioncount;
	for (int i = 0; i < regioncount; i++)
	{
		regions[i].w = regionwidth;
		regions[i].h = SIMULATIONH;
		regions[i].x = regionwidth * i;
		regions[i].y = 0;

		if (i == regioncount - 1)
		{
			if ((regions[i].w + regions[i].x) != SIMULATIONW)
			{
				regions[i].w += SIMULATIONW - (regions[i].w + regions[i].x);
			}
		}

		region_groups[i % region_group_count][i / region_group_count] = regions[i];
	}

	std::cout << threadcount << " threads, " << region_group_count << " region groups, " << regioncount << " regions" << std::endl;
}

void reinit_simulation(int threadcount_, int groupcount_, atom *parts)
{
	exiting = true;
	for (int i = 0; i < threadcount; i++)
	{
		locks[i]->clear();
		threads[i].join();
	}
	exiting = false;

	for (int i = 0; i < threadcount; i++)
	{
		delete locks[i];
	}
	delete[] locks;
	delete[] threads;

	for (int i = 0; i < region_group_count; i++)
		delete[] region_groups[i];
	delete[] region_groups;
	delete[] regions;

	init_simulation(threadcount_, groupcount_, parts);
}

void simulate()
{
	last_partcount = 0;

	for (int j = 0; j < region_group_count; j++)
	{
		active_regions = region_groups[j];

		runs = 0;

		for (int i = 0; i < threadcount; i++)
		{
			locks[i]->clear();
		}

		while (runs < threadcount)
			;
	}

	mutex = !mutex;
}

void add_parts(atom *parts, int origin_x, int origin_y, uint8_t type)
{
	int radius = 10;
	for (int y = origin_y - radius; y < origin_y + radius; y++)
	{
		if (y < 0 || y >= SIMULATIONH)
		{
			continue;
		}
		for (int x = origin_x - radius; x < origin_x + radius; x++)
		{
			if (x < 0 || x >= SIMULATIONW)
			{
				continue;
			}
			parts[PART(x, y)].type = type;
			parts[PART(x, y)].vx = 0;
			parts[PART(x, y)].vy = 0;
			parts[PART(x, y)].x = x;
			parts[PART(x, y)].y = y;
			if (type == TYPE_PARTICLE)
			{
				parts[PART(x, y)].vx = randfd() * 5.0f;
				parts[PART(x, y)].vy = randfd() * 5.0f;
			}
		}
	}
}

void fill_screen(atom *parts, uint8_t type)
{
	for (int x = 0; x < SIMULATIONW / 2; x++)
	{
		for (int y = 0; y < SIMULATIONH; y++)
		{
			parts[PART(x, y)].type = type;
			parts[PART(x, y)].vy = 0;
			parts[PART(x, y)].x = x;
			parts[PART(x, y)].y = y;
			if (type == TYPE_PARTICLE)
			{
				parts[PART(x, y)].vx = randfd() * 5.0f;
				parts[PART(x, y)].vy = randfd() * 5.0f;
			}
		}
	}
}

int main(int argc, char const *argv[])
{
	int num_threads = 4;
	int num_groups = 2;

	if(argc > 2)
	{
		try {
			num_threads = std::stoi(argv[1]);
			num_groups = std::stoi(argv[2]);
		} catch (std::exception &e) {
			std::cout << argv[0] << " num_threads num_groups" << std::endl;
			return -1;
		}
	}

	atom *parts = new atom[SIMULATIONW * SIMULATIONH];
	std::fill(parts, parts + (SIMULATIONW * SIMULATIONH), atom());

	init_simulation(num_threads, num_groups, parts);

	uint8_t particle_type = TYPE_POWDER;
	fill_screen(parts, particle_type);

	long sim_steps = 0;

	auto start = std::chrono::high_resolution_clock::now();

	while (true)
	{
		simulate();
		sim_steps++;

		if (sim_steps == 500)
		{
			break;
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	std::cout << sim_steps << " frames in " << total_time.count() << "ms, " << total_time.count() / (double)sim_steps << "ms per frame" << std::endl;

	return 0;
}