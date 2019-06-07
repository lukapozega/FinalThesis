#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream> 
#include <sstream>
#include <tuple>
#include <algorithm>
#include <unordered_map>

#include "Config.h"

#include "bioparser/bioparser.hpp"
#include "PAFObject.cpp"
#include "FASTAQObject.cpp"
#include "repeats_parser.h"

struct option options[] = {
	{"version", no_argument, 0, 'v'},
	{"help", no_argument, 0, 'h'},
	{0, 0, 0, 0}
};

std::unordered_map<char, char> complement_map = {
			{'A', 'T'},
			{'T', 'A'},
			{'G', 'C'},
			{'C', 'G'},
	};

class Vertex {

public:
	PAFObject read;
	std::vector<Vertex*> vertices;
	Vertex* parent = NULL;

	Vertex(PAFObject read) : read(read) {
		this->read = read;
	}
};

void add_breakpoints(std::vector<std::unique_ptr<PAFObject>> &paf_objects, std::vector<std::tuple<std::string, int, int>> &repeats) {
	std::vector<std::unique_ptr<PAFObject>>::iterator it = paf_objects.begin();
	std::vector<std::unique_ptr<PAFObject>>::iterator to_remove;
	while (it != paf_objects.end()) {
		for (auto const& r : repeats) {
			if ((*it)->t_begin > std::get<1>(r) && (*it)->t_end < std::get<2>(r)) {
				to_remove = it++;
				(*to_remove) = NULL;
				if (it == paf_objects.end()) {
					paf_objects.erase(std::remove_if(paf_objects.begin(), paf_objects.end(), [](std::unique_ptr<PAFObject> &p) {return p == NULL;}), paf_objects.end());
					return;
				}
			}
		}
		it++;
	}
	paf_objects.erase(std::remove_if(paf_objects.begin(), paf_objects.end(), [](std::unique_ptr<PAFObject> &p) {return p == NULL;}), paf_objects.end());
}

std::unordered_map<std::string, std::vector<Vertex*>> create_graph(std::vector<Vertex> &vertices, std::vector<std::unique_ptr<PAFObject>> &paf_objects) {
	for (auto const& paf: paf_objects) {
		vertices.emplace_back(Vertex(*paf));
	}
	std::vector<Vertex>::iterator it = vertices.begin();
	std::vector<Vertex>::iterator next;
	std::vector<Vertex*> heads;
	std::unordered_map<std::string, std::vector<Vertex*>> re;
	heads.emplace_back(&(*it));
	while(it != --vertices.end()) {
		next = std::next(it);
		if ((*next).read.t_name != (*it).read.t_name) {
			re[(*it).read.t_name] = heads;
			heads.clear();
			heads.emplace_back(&(*next));
			it++;
			continue;
		}
		while ((*next).read.t_begin <= (*it).read.t_end && (*next).read.t_name == (*it).read.t_name) {
			it->vertices.emplace_back(&(*next));
			next++;
			if(next == vertices.end()) break;
		}
		if (next == std::next(it)){
		 	heads.emplace_back(&(*next));
		}
		it++;
	}
	re[(*it).read.t_name] = heads;
	heads.clear();
	FILE* file = fopen ("network.gfa","w");
	for (auto const& vertex : vertices) {
		fprintf (file, "S\t%s\t%c\tLN:i:%u\n", vertex.read.q_name.c_str(), '*', vertex.read.q_length);
		for (auto const& next: vertex.vertices) {
			fprintf (file, "L\t%s\t%c\t%s\t%c\t%c\n", vertex.read.q_name.c_str(), vertex.read.orientation, next->read.q_name.c_str(), next->read.orientation, '*');
		}
	}
	fclose (file);
	return re;
}

std::vector<Vertex*> DepthFirstSearch(std::unordered_map<std::string, std::vector<Vertex*>> heads) {
	Vertex* max;
	Vertex* head;
	Vertex* longest;
	std::vector<Vertex*> ends;
	int begin;
	for (auto const& seq: heads) {
		int length=0;
		for (int i = 0; i < seq.second.size(); i++) {
			head = seq.second[i];
			begin = head->read.t_begin;
			while(!head->vertices.empty()) {
				max = head->vertices[0];
				for (auto const& vertex : head->vertices) {
					if (vertex->read.t_begin > max->read.t_begin) max = vertex;
				}
				max->parent = head;
				head = max;
			}
			if (length < head->read.t_begin - begin) {
				length = head->read.t_begin - begin;
				longest = head;
			}
		}
		ends.emplace_back(longest);
	}
	return ends;
}

bool paf_unique(const std::unique_ptr<PAFObject>& a, const std::unique_ptr<PAFObject>& b) {
	return a->q_name == b->q_name;
}

void clear_contained_reads(std::vector<std::unique_ptr<PAFObject>> &paf_objects) {
	auto long_cmp = [](const std::unique_ptr<PAFObject>& a, const std::unique_ptr<PAFObject>& b) {
        if (a->q_name == b->q_name) {
            return (a->t_end - a->t_begin > b->t_end - b->t_begin);
        }
        return (a->q_name > b->q_name);
    };
    std::sort(paf_objects.begin(), paf_objects.end(), long_cmp);
	std::vector<std::unique_ptr<PAFObject>>::iterator it = paf_objects.begin();
	it = std::unique (paf_objects.begin(), paf_objects.end(), paf_unique);
	paf_objects.resize(std::distance(paf_objects.begin(), it));
	paf_objects.erase(std::remove_if(paf_objects.begin(), paf_objects.end(), [](std::unique_ptr<PAFObject> &p){return p->q_end - p->q_begin < 0.9 * p->q_length;}), paf_objects.end());
	auto paf_cmp = [](const std::unique_ptr<PAFObject>& a, const std::unique_ptr<PAFObject>& b) {
		if (a->t_name == b->t_name) { 
	 		if (a->t_begin == b->t_begin) {
	        	return (a->t_end > b->t_end);
	        }
	        return (a->t_begin < b->t_begin);
	    }
	    return (a->t_name > b->t_name);
	};
	std::sort(paf_objects.begin(), paf_objects.end(), paf_cmp);
	it = paf_objects.begin();
	std::vector<std::unique_ptr<PAFObject>>::iterator next;
	std::vector<std::unique_ptr<PAFObject>>::iterator to_remove;
	while (it != --paf_objects.end()) {
        do {
        	next = std::next(it);
        } while ((*next) == NULL && (*next)->t_name == (*it)->t_name);
        while ((*next)->t_end <= (*it)->t_end && (*next)->t_name == (*it)->t_name) {
            to_remove = next;
            do {
                next++;
            } while ((*next) == NULL && next != paf_objects.end());
            (*to_remove) = NULL;
            if(next == paf_objects.end()) {
                paf_objects.erase(std::remove_if(paf_objects.begin(), paf_objects.end(), [](std::unique_ptr<PAFObject> &p) {return p == NULL;}), paf_objects.end());
                return;
            }
        }
        it = next;
    }
    paf_objects.erase(std::remove_if(paf_objects.begin(), paf_objects.end(), [](std::unique_ptr<PAFObject> &p) {return p == NULL;}), paf_objects.end());
}

bool file_format(const std::string &str, const std::string &suffix) {
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string complement(std::string const &sequence){
	std::string comp = "";
	for (char const & c : sequence) {
		comp = complement_map.at(c) + comp; 
	}
	return comp;
}

void statistics(std::vector<Vertex*> ends, std::vector<std::unique_ptr<FASTAQEntity>> & ref_objects, std::string const &file_path) {
	std::vector<std::unique_ptr<FASTAQEntity>> reads;
	if (file_format(file_path, ".fastq") || file_format(file_path, ".fastq.gz")) {
		auto reads_parser = bioparser::createParser<bioparser::FastqParser, FASTAQEntity>(file_path);
		reads_parser->parse(reads, -1);
	} else {
		auto reads_parser = bioparser::createParser<bioparser::FastaParser, FASTAQEntity>(file_path);
		reads_parser->parse(reads, -1);
	}
	std::string part;
	printf("Writing to files...\n\n");
	FILE* file = fopen ("used_reads.txt","w");
	FILE* genome = fopen ("genome.fasta","w");
	int startIndex, endIndex, last_index, overlap, used_reads;
	std::string seq;
	Vertex* curr;
	for (auto const& end : ends) {
		startIndex = end->read.q_begin;
		endIndex = end->read.q_end;
		last_index = end->read.t_end;
		overlap = 0;
		used_reads = 0;
		curr = end;
		fprintf (genome, ">%s\n", end->read.t_name.c_str());
		while (curr->parent != NULL) {
			used_reads++;
			for (auto const& r : reads) {
				if (r->name == curr->read.q_name) {
					if (curr->read.orientation == '+') {
						seq = complement(r->sequence);
					} else {
						seq = r->sequence;
					}
					part.assign(seq, startIndex + overlap, endIndex - (startIndex + overlap));
					break;
				}
			}
			fprintf (genome, "%s", part.c_str());
			endIndex = curr->read.q_begin;
			fprintf (file, "%s %s %d %d\n", curr->read.q_name.c_str(), curr->read.t_name.c_str(), curr->read.t_begin, curr->read.t_end);
			overlap = curr->read.t_begin;
			curr = curr->parent;
			overlap = curr->read.t_end - overlap;
			startIndex = curr->read.q_begin;
		}
		fprintf(genome, "\n");
		printf("%s coverage: %f%%\n", end->read.t_name.c_str(), (last_index-curr->read.t_begin) / (float) curr->read.t_length * 100);
		printf("Number of used reads: %d\n",used_reads);
	}
	fclose (file);
	fclose (genome);
}

void help() {
	printf("Program accepts three arguments and prints assemblying statistics.\n");
	printf("Usage: assembly <alignment file > <reference file> <repeats file> <reads file>\n");
	printf("Arguments should be in PAF, fasta rpt and fasta/fastq file format\n");
}

void version() {
	printf("Version %d.%d\n", assembly_VERSION_MAJOR, assembly_VERSION_MINOR);
}

int main(int argc, char** argv) {

	char optchr;
	int option_index = 0;
	while((optchr = getopt_long(argc, argv, "hv", options, &option_index)) != -1) {
		switch (optchr) {
			case 'h': {
				help();
				return 1;
			}
			case 'v': {
				version();
				return 1;
			}
			default: {
				fprintf(stderr, "Unknown option -%c\n", optchr);
				return 1;
			}
		}
	}

	if(argc - optind < 4) {
		fprintf(stderr, "Program requires three arguments.\n");
		fprintf(stderr, "Use \"-h\" or \"--help\" for more information.\n");
		return 1;
	}

	std::vector<std::unique_ptr<PAFObject>> paf_objects;
	auto paf_parser = bioparser::createParser<bioparser::PafParser, PAFObject>(argv[optind]);
	paf_parser->parse(paf_objects, -1);
	clear_contained_reads(paf_objects);

	std::vector<std::unique_ptr<FASTAQEntity>> ref_objects;
	auto fasta_parser = bioparser::createParser<bioparser::FastaParser, FASTAQEntity>(argv[optind + 1]);
	fasta_parser->parse(ref_objects, -1);

	std::vector<std::tuple<std::string, int, int>> repeats;
	if (!repeats_parser::parse(repeats, argv[3])) {
		fprintf(stderr, "Error reading file %s\n", argv[optind + 2]);
		return 1;
	}

	repeats_parser::remove_covered(repeats, paf_objects);

	repeats_parser::check_repeats(repeats, ref_objects);

	add_breakpoints(paf_objects, repeats);

	std::vector<Vertex> vertices;
	std::unordered_map<std::string, std::vector<Vertex*>> heads = create_graph(vertices, paf_objects);

	std::vector<Vertex*> ends = DepthFirstSearch(heads);

	statistics(ends, ref_objects, argv[optind + 3]);

	return 0;

}


