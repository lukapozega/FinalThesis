#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream> 
#include <sstream>
#include <tuple>
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
		if ((*next).read.target_name != (*it).read.target_name) {
			re[(*it).read.target_name] = heads;
			heads.clear();
			heads.emplace_back(&(*next));
			it++;
			continue;
		}
		while ((*next).read.t_begin <= (*it).read.t_end && (*next).read.target_name == (*it).read.target_name) {
			it->vertices.emplace_back(&(*next));
			next++;
			if(next == vertices.end()) break;
		}
		if (next == std::next(it) && (*next).read.target_name == (*std::next(it)).read.target_name) heads.emplace_back(&(*next));
		it++;
	}
	re[(*it).read.target_name] = heads;
	heads.clear();
	// FILE* file = fopen ("network.txt","w");
	// for (auto const& vertex : vertices) {
	// 	fprintf (file, "%d %d", vertex.read.t_begin, vertex.read.t_end);
	// 	for (auto const& next: vertex.vertices) {
	// 		fprintf (file, " %d", next->read.t_begin);
	// 	}
	// 	fprintf (file, "\n");
	// }
	// fclose (file);
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
	return a->query_name == b->query_name;
}

void clear_contained_reads(std::vector<std::unique_ptr<PAFObject>> &paf_objects) {
	std::vector<std::unique_ptr<PAFObject>>::iterator it = paf_objects.begin();
	it = std::unique (paf_objects.begin(), paf_objects.end(), paf_unique);
	paf_objects.resize(std::distance(paf_objects.begin(), it));
	paf_objects.erase(std::remove_if(paf_objects.begin(), paf_objects.end(), [](std::unique_ptr<PAFObject> &p){return p->q_end - p->q_begin < 0.9 * p->q_length;}), paf_objects.end());
	auto paf_cmp = [](const std::unique_ptr<PAFObject>& a, const std::unique_ptr<PAFObject>& b) {
		if (a->target_name == b->target_name) { 
	 		if (a->t_begin == b->t_begin) {
	        	return (a->t_end > b->t_end);
	        }
	        return (a->t_begin < b->t_begin);
	    }
	    return (a->target_name > b->target_name);
	};
	std::sort(paf_objects.begin(), paf_objects.end(), paf_cmp);
	it = paf_objects.begin();
	int s, e;
	std::string n;
	while (it != --paf_objects.end()) {
		s = (*it)->t_begin;
		e = (*it)->t_end;
		n = (*it)->target_name;
		paf_objects.erase(std::remove_if(it+1, paf_objects.end(), [&s, &e, &n](std::unique_ptr<PAFObject> &p){return p->t_begin >= s && p->t_end <= e && p->target_name == n;}), paf_objects.end());
		if (it != --paf_objects.end()) {
			it++;
		}
	}
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
	if (file_format(file_path, ".fastq")) {
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
		fprintf (genome, ">%s\n", end->read.target_name.c_str());
		while (curr->parent != NULL) {
			used_reads++;
			for (auto const& r : reads) {
				if (r->name == curr->read.query_name) {
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
			fprintf (file, "%s %s %d %d\n", curr->read.query_name.c_str(), curr->read.target_name.c_str(), curr->read.t_begin, curr->read.t_end);
			overlap = curr->read.t_begin;
			curr = curr->parent;
			overlap = curr->read.t_end - overlap;
			startIndex = curr->read.q_begin;
		}
		fprintf(genome, "\n");
		printf("%s coverage: %f%%\n", end->read.target_name.c_str(), (last_index-curr->read.t_begin) / (float) curr->read.t_length * 100);
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

	std::vector<Vertex> vertices;
	std::unordered_map<std::string, std::vector<Vertex*>> heads = create_graph(vertices, paf_objects);

	std::vector<Vertex*> ends = DepthFirstSearch(heads);

	statistics(ends, ref_objects, argv[optind + 3]);

	return 0;

}


