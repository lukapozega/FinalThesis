namespace repeats_parser {
			   
	bool parse(std::vector<std::tuple<std::string, int, int>> &repeats, const std::string& path);

	void remove_covered(std::vector<std::tuple<std::string, int, int>> &repeats, std::vector<std::unique_ptr<PAFObject>> &paf_objects);

	std::tuple<int, int, int, int> check_repeats(std::vector<std::tuple<std::string, int, int>> &repeats, const std::string& reference);

}