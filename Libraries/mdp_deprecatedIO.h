/////////////////////////////////////////////////////////////////
/// @file mdp_deprecatedIO.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Old functions for file IO now deprecated
///
/// Distributed under GPL2 License
/// Read attached license in file mdp_license.txt
/// This file cannot be distributed without file mdp_license.txt
//////////////////////////////////////////////////////////////////
#ifndef MDP_DEPRECATEDIO_
#define MDP_DEPRECATEDIO_

namespace MDP
{
#ifdef INCLUDE_DEPRECATED_IO
	template <class T>
	void mdp_field<T>::load(const char filename[],
							int processIO,
							mdp_int max_buffer_size,
							const char *header,
							mdp_int header_size,
							mdp_int (*sort_x)(mdp_lattice &, mdp_int),
							int auto_switch_endianess)
	{
		mdp_int idx_gl, nvol_gl = lattice().global_volume();
		double mytime = mdp.time();
		int try_switch_endianess = false;
		if (isSubProcess(processIO))
		{
			auto buffer_size = std::make_unique<mdp_int[]>(Nproc);
			mdp_array<T, 3> large_buffer(Nproc, max_buffer_size, m_field_components);
			auto short_buffer = std::make_unique<T[]>(m_field_components);
			int process;
			mdp_request request;

			for (process = 0; process < Nproc; process++)
				buffer_size[process] = 0;
			printf("Loading file %s from process %i (buffer = %li sites)\n",
				   filename, processIO, max_buffer_size);
			fflush(stdout);
			FILE *fp = fopen(filename, "rb");
			if (fp == nullptr)
				error("Unable to open file");

			// #ifdef INSERT_HEADER
			if (strcmp(header, "NATIVE") == 0)
			{
				error("NATIVE HEADER IN DEPRECATED FUNCTION NOT SUPPORTED ANY MORE");
			}
			else if (header != nullptr && strcmp(header, "NOHEADER") != 0)
			{
				if (fread(header, sizeof(char), header_size, fp) !=
					header_size)
					error("Unable to read file header");
			}
			// #endif
			for (idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
			{
				process = where_global(idx_gl);
				if (process != NOWHERE)
				{
					if (sort_x != 0)
						if (fseek(fp, sort_x(lattice(), idx_gl) * sizeof(T) * m_field_components + header_size, SEEK_SET) < 0)
							error("unexpected end of file");
					if ((fread(short_buffer.get(), sizeof(T), m_field_components, fp) -
						 m_field_components) != 0)
						error("unexpected end of file");
				}
				if ((process != NOWHERE) && (process != processIO))
				{
					for (mdp_uint k = 0; k < m_field_components; k++)
						large_buffer(process, buffer_size[process], k) = short_buffer[k];
					buffer_size[process]++;
					if (buffer_size[process] == max_buffer_size)
					{
						mdp.put(&(large_buffer(process, 0, 0)),
								max_buffer_size * m_field_components, process, request);
						mdp.wait(request);
						buffer_size[process] = 0;
					}
					if (idx_gl == nvol_gl - 1)
						for (process = 0; process < Nproc; process++)
							if ((process != ME) &&
								(buffer_size[process] != max_buffer_size) &&
								(buffer_size[process] > 0))
							{
								mdp.put(&(large_buffer(process, 0, 0)),
										buffer_size[process] * m_field_components,
										process, request);
								mdp.wait(request);
							}
				}
				if (process == processIO)
				{
					for (mdp_uint k = 0; k < m_field_components; k++)
						*(m_data.get() + lattice().local(idx_gl) * m_field_components + k) = short_buffer[k];
				}
			}
			fclose(fp);
		}
		else
		{
			mdp_int buffer_size = 0;
			auto local_index = std::make_unique<mdp_int[]>(max_buffer_size);
			mdp_array<T, 2> local_buffer(max_buffer_size, m_field_components);
			for (idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
			{
				int process = where_global(idx_gl);
				if (process == ME)
				{
					local_index[buffer_size] = lattice().local(idx_gl);
					buffer_size++;
				}

				if ((buffer_size == max_buffer_size) || ((idx_gl == nvol_gl - 1) && (buffer_size > 0)))
				{
					mdp.get(&(local_buffer(0, 0)), buffer_size * m_field_components, processIO);
					for (mdp_int idx = 0; idx < buffer_size; idx++)
						for (mdp_uint k = 0; k < m_field_components; k++)
							m_data[local_index[idx] * m_field_components + k] = local_buffer(idx, k);

					buffer_size = 0;
				}
			}
		}

		update();

		if (isMainProcess() && !mdp_shutup)
		{
			printf("... Loading time: %f (sec)\n", mdp.time() - mytime);
			fflush(stdout);
		}

		if (try_switch_endianess && auto_switch_endianess)
		{
			if (isMainProcess())
				printf("Trying to switch endianess.\n");
			switch_endianess_4bytes();
		}
	}

	template <class T>
	void mdp_field<T>::save(const char filename[],
							int processIO,
							mdp_int max_buffer_size,
							const char *header,
							mdp_int header_size,
							mdp_int (*sort_x)(mdp_lattice &, mdp_int),
							const char *mode)
	{
		mdp_int nvol_gl = lattice().global_volume();
		double mytime = mdp.time();

		if (isSubProcess(processIO))
		{
			auto buffer_size = std::make_unique<mdp_int[]>(Nproc);
			auto buffer_ptr = std::make_unique<mdp_int[]>(Nproc);
			mdp_array<T, 3> large_buffer(Nproc, max_buffer_size, m_field_components);
			auto short_buffer = std::make_unique<T[]>(m_field_components);

			for (int process = 0; process < Nproc; process++)
				buffer_ptr[process] = 0;
			printf("Saving file %s from process %i (buffer = %li sites)\n",
				   filename, processIO, max_buffer_size);
			fflush(stdout);
			FILE *fp = fopen(filename, mode);
			if (fp == nullptr)
				error("Unable to open file");

			// Write optional header
			if (strcmp(header, "NATIVE") == 0)
			{
				error("NATIVE HEADER IN DEPRECATED FUNCTION NOT SUPPORTED ANY MORE");
			}
			else if (header != nullptr && strcmp(header, "NOHEADER") != 0)
			{
				if (fwrite(header, sizeof(char), header_size, fp) !=
					header_size)
					error("Unable to write file header");
			}

			// Loop over global sites
			for (mdp_int idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
			{
				int process = where_global(idx_gl);

				// Receive data from other processes
				if ((process != NOWHERE) && (process != processIO))
				{
					if (buffer_ptr[process] == 0)
					{
						mdp.get(buffer_size[process], process);
						mdp.get(&(large_buffer(process, 0, 0)),
								buffer_size[process] * m_field_components, process);
					}

					for (mdp_uint k = 0; k < m_field_components; k++)
						short_buffer[k] = large_buffer(process, buffer_ptr[process], k);

					buffer_ptr[process]++;
					if (buffer_ptr[process] == buffer_size[process])
						buffer_ptr[process] = 0;
				}

				// Local data
				if (process == processIO)
				{
					mdp_int local_idx = lattice().local(idx_gl);
					for (mdp_uint k = 0; k < m_field_components; k++)
						short_buffer[k] = m_data[local_idx * m_field_components + k];
				}

				// Write to file
				if (process != NOWHERE)
				{
					if (sort_x != 0)
						if (fseek(fp, sort_x(lattice(), idx_gl) * sizeof(T) * m_field_components + header_size, SEEK_SET) < 0)
							error("unexpected end of file");
					if ((fwrite(short_buffer.get(), sizeof(T), m_field_components, fp) -
						 m_field_components) != 0)
						error("I cannot write on the file. I am confused !?!?");
				}
			}

			if (strcmp(header, "NATIVE") == 0)
				fprintf(fp, "\n\n [ MDP Standard File Format ]\n");

			fclose(fp);
		}
		else // Main process
		{
			mdp_int buffer_size = 0;
			auto local_index = std::make_unique<mdp_int[]>(max_buffer_size);
			mdp_array<T, 2> local_buffer(max_buffer_size, m_field_components);
			mdp_request request;

			for (mdp_int idx_gl = 0; idx_gl < nvol_gl; idx_gl++)
			{
				int process = where_global(idx_gl);
				if (process == ME)
				{
					local_index[buffer_size] = lattice().local(idx_gl);
					buffer_size++;
				}

				if ((buffer_size == max_buffer_size) || ((idx_gl == nvol_gl - 1) && (buffer_size > 0)))
				{
					for (mdp_int idx = 0; idx < buffer_size; idx++)
						for (mdp_uint k = 0; k < m_field_components; k++)
							local_buffer(idx, k) = m_data[local_index[idx] * m_field_components + k];

					mdp.put(buffer_size, processIO, request);
					mdp.wait(request);
					mdp.put(&(local_buffer(0, 0)), buffer_size * m_field_components,
							processIO, request);
					mdp.wait(request);
					buffer_size = 0;
				}
			}
		}

		if (isMainProcess() && !mdp_shutup)
		{
			printf("... Saving time: %f (sec)\n", mdp.time() - mytime);
			fflush(stdout);
		}
	}
#endif
} // namespace MDP

#endif /* MDP_DEPRECATEDIO_ */
