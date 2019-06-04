class PAFObject {

public:

	std::string q_name;
	std::uint32_t q_length;
	std::uint32_t q_begin;
	std::uint32_t q_end;
	char orientation;
	std::string t_name;
	std::uint32_t t_length;
	std::uint32_t t_begin;
	std::uint32_t t_end;

	PAFObject(
		const char* q_name, std::uint32_t q_name_length,
	    std::uint32_t q_length,
	    std::uint32_t q_begin,
	    std::uint32_t q_end,
	    char orientation,
	    const char* t_name, std::uint32_t t_name_length,
	    std::uint32_t t_length,
	    std::uint32_t t_begin,
	    std::uint32_t t_end,
	    std::uint32_t matching_bases,
	    std::uint32_t overlap_length,
	    std::uint32_t mapping_quality) {
			(this->q_name).assign(q_name, q_name_length);
			this->q_length = q_length;
			this->q_begin = q_begin;
			this->q_end = q_end;
			this->orientation = orientation;
			(this->t_name).assign(t_name, t_name_length);
			this->t_length = t_length;
			this->t_begin = t_begin;
			this->t_end = t_end;
    }

};