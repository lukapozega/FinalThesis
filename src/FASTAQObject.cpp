class FASTAQEntity {
	
public:
	std::string name;
	std::string sequence;
	
	FASTAQEntity(
		const char *name, uint32_t name_length,
		const char *sequence, uint32_t sequence_length) {
		(this->name).assign(name, name_length);
		(this->sequence).assign(sequence, sequence_length);
	}

	FASTAQEntity (
		const char* name, std::uint32_t name_length,
        const char* sequence, std::uint32_t sequence_length,
        const char* quality, std::uint32_t quality_length) {
		(this->name).assign(name, name_length);
		(this->sequence).assign(sequence, sequence_length);
	}
};