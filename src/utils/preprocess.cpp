#include "preprocess.h"

void preprocessInputFasta(const std::string input_path, const std::string workdir) {
    spdlog::info("Preprocessing input FASTA file: {}", input_path);
    spdlog::info("Working directory: {}", workdir);

    // 确保工作目录下有data文件夹
    file_io::FilePath data_dir = file_io::FilePath(workdir) / WORKDIR_DATA;
    file_io::ensureDirectoryExists(data_dir);

    // 在data文件夹下创建raw_data文件夹
    file_io::FilePath raw_data_dir = data_dir / DATA_RAW;
    file_io::ensureDirectoryExists(raw_data_dir);

    // 在data文件夹下创建clean_data文件夹
    file_io::FilePath clean_data_dir = data_dir / DATA_CLEAN;
    file_io::ensureDirectoryExists(clean_data_dir);

    // 将inputdata复制到raw_data文件夹
    file_io::FilePath input_file = file_io::FilePath(input_path);
    file_io::FilePath dest_file = raw_data_dir / input_file.filename();

    if (file_io::isUrl(input_file))
    {
        file_io::downloadFile(input_file, dest_file);
    }else
    {
        file_io::copyFile(input_file, dest_file);
    }

    spdlog::info("Preprocessing completed.");
}