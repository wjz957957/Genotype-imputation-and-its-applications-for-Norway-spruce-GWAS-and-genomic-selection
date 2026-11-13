import gzip
import numpy as np
import subprocess


# 定义阈值范围列表，每个范围是一个元组(min_threshold, max_threshold)
threshold_ranges = [
    (0.00, 0.05),
    (0.05, 0.10),
    (0.10, 0.15),
    (0.15, 0.20),
    (0.20, 0.25),
    (0.25, 0.30),
    (0.30, 0.35),
    (0.35, 0.40),
    (0.40, 0.45),
    (0.45, 0.50),
    (0.50, 0.55), 
    (0.55, 0.60), 
    (0.60, 0.65), 
    (0.65, 0.70), 
    (0.70, 0.75), 
    (0.75, 0.80), 
    (0.80, 0.85), 
    (0.85, 0.90), 
    (0.90, 0.95), 
    (0.95, 1.01)
]
# 遍历每个阈值范围
for min_threshold, max_threshold in threshold_ranges:
    # 构建bcftools view命令
    bcftools_view_command = [
        'bcftools', 'view', '-i',
        f"INFO/INFO>={min_threshold}&INFO/INFO<{max_threshold}",
        'impute5_output_map.vcf.gz'
    ]
    
    # 定义bgzip命令
    bgzip_command = ['bgzip', '-c']
    
    # 将bcftools的输出重定向到bgzip
    # 使用subprocess.Popen()创建bcftools进程，并将输出设置为管道
    bcftools_process = subprocess.Popen(bcftools_view_command, stdout=subprocess.PIPE, text=True)
    
    # 将bcftools的输出作为bgzip的输入
    with open('impute5_output_INFO-0.1-concordance.vcf.gz', 'w') as output_file:
        bgzip_process = subprocess.Popen(bgzip_command, stdin=bcftools_process.stdout, stdout=output_file)
    
    # 等待bcftools进程完成
    bcftools_process.wait()
    
    # 确保bcftools进程正常结束
    if bcftools_process.returncode != 0:
        print("bcftools encountered an error")
    
    # 等待bgzip进程完成
    bgzip_process.wait()
    
    # 确保bgzip进程正常结束
    if bgzip_process.returncode != 0:
        print("bgzip encountered an error")
    	
    ##提取筛选后的染色体位置信息	
    
    def extract_columns_vcf(input_vcf_gz_path, output_file_path):
        with gzip.open(input_vcf_gz_path, 'rt') as infile, open(output_file_path, 'w') as outfile:
            for line in infile:
                # 跳过以'#'开头的注释行
                if line.startswith('#'):
                    continue
                # 分割行并提取第一列和第二列
                columns = line.strip().split('\t')
                if len(columns) >= 2:
                    # 写入到输出文件
                    outfile.write(f"{columns[0]}\t{columns[1]}\n")
    
    # 使用示例
    input_vcf_gz_path = 'impute5_output_INFO-0.1-concordance.vcf.gz'  # 你的.vcf.gz文件路径
    output_file_path = 'extracted_chr-pos-concordance.txt'    # 你希望输出的文件路径
    extract_columns_vcf(input_vcf_gz_path, output_file_path)
    
    ##根据提出的染色体位置信息，将重测序里相同的位点提取出来
    
    # 定义bcftools view命令和参数
    bcftools_command = [
        'bcftools', 'view', '--regions-file', 'extracted_chr-pos-concordance.txt',
        'Chr11_map_phased_50_rename.vcf.gz'
    ]
    
    # 定义bgzip命令
    bgzip_command = ['bgzip', '-c']
    
    # 创建bcftools进程，并将输出重定向到bgzip的输入
    with open('Chr11_map_phased_50_rename_INFO-0.1-concordance.vcf.gz', 'wb') as output_file:
        bcftools_process = subprocess.Popen(bcftools_command, stdout=subprocess.PIPE, text=False)
        bgzip_process = subprocess.Popen(bgzip_command, stdin=bcftools_process.stdout, stdout=output_file)
    
    # 等待bcftools进程完成
    bcftools_process.wait()
    
    # 确保bcftools进程正常结束
    if bcftools_process.returncode != 0:
        print("bcftools encountered an error")
    
    # 等待bgzip进程完成
    bgzip_process.wait()
    
    # 确保bgzip进程正常结束
    if bgzip_process.returncode != 0:
        print("bgzip encountered an error")




    #######################提取参考文件GT字段
    # 指定VCF文件路径
    vcf_filename = 'Chr11_map_phased_50_rename_INFO-0.1-concordance.vcf.gz'
    
    # 定义基因型到数组的映射
    genotype_to_array = {
        '0|0': '0',
        '1|0': '1',
        '0|1': '1',
        '1|1': '2'
    }
    
    # 用于存储转换后的数据列表
    converted_data_list = []
    
    # 定义数据转换函数
    def convert_genotypes(genotypes_str):
        # 根据提供的映射转换基因型字符串，并用制表符分隔结果
        return '\t'.join(genotype_to_array[gt] for gt in genotypes_str.split('\t') if gt in genotype_to_array)
    
    # 使用 gzip 打开文件并读取内容
    with gzip.open(vcf_filename, 'rt') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue  # 跳过注释行
            columns = line.strip().split('\t')
            if len(columns) > 9:
                # 提取基因型数据，假设每行有50个基因型数据
                genotypes_str = '\t'.join(columns[9:])  # 根据VCF格式提取基因型字段
                # 转换基因型数据
                converted_data = convert_genotypes(genotypes_str)
                if converted_data:  # 确保转换结果不为空
                    converted_data_list.append(converted_data)
    
    # 将转换后的数据写入新文件
    output_filename = 'converted_genotypes-concordance.txt'
    with open(output_filename, 'w') as output_file:
        for converted_data in converted_data_list:  # 确保我们使用的是正确的变量名
            output_file.write(converted_data + '\n')
    
    print("Conversion complete. Results are saved in", output_filename)
    
    
    
    ########################################提取填充文件GT字段
    # 指定VCF文件路径
    vcf_filename = 'impute5_output_INFO-0.1-concordance.vcf.gz'  # 替换为您的文件名
    
    # 用于存储提取的数据
    extracted_data = []
    
    # 使用 gzip 打开文件并读取内容
    with gzip.open(vcf_filename, 'rt') as vcf_file:
        for line in vcf_file:
            # 跳过注释行
            if line.startswith('#'):
                continue
            # 拆分行数据
            columns = line.strip().split('\t')
            # 检查列数是否足够
            if len(columns) > 9:
                # 提取第十列及以后的数据
                extracted_columns = '\t'.join(columns[9:])
                # 存储提取的数据
                extracted_data.append(extracted_columns)
    
    # 将提取的数据写入到新文件
    output_filename = 'extracted_columns-imputed-concordance.txt'  # 您可以自定义输出文件名
    with open(output_filename, 'w') as output_file:
        for data in extracted_data:
            output_file.write(data + '\n')
    
    # 打印完成消息
    print("Extraction complete. Results are saved in", output_filename)
    
    # 指定输入和输出文件路径
    input_filename = 'extracted_columns-imputed-concordance.txt'
    output_filename = 'extracted_values-imputed-concordance.txt'
    
    # 读取输入文件并写入输出文件
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        # 逐行读取文件
        for line_number, line in enumerate(infile):
            
            # 处理第二行及以后的行
            fields = line.strip().split('\t')  # 根据制表符分割行
            extracted_values = [field.split(':')[2] for field in fields if ':' in field]  # 提取冒号后的值
            extracted_line = '\t'.join(extracted_values)  # 将提取的值重新连接成一行
            outfile.write(extracted_line + '\n')  # 写入结果到输出文件
    
    # 打印完成消息
    print("Extraction complete. Results are saved in", output_filename)
    
    
    
    ##########合并两个文件
    def merge_files(file1_path, file2_path, output_path):
        try:
            # 打开两个文件并读取内容
            with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
                file1_lines = file1.readlines()
                file2_lines = file2.readlines()
    
            # 确保两个文件的行数相同
            if len(file1_lines) != len(file2_lines):
                raise ValueError("两个文件的行数必须相同")
    
            # 合并单元格并连接数据
            merged_lines = []
            for line1, line2 in zip(file1_lines, file2_lines):
                # 移除行尾的换行符并按制表符分割
                cells1 = line1.strip().split('\t')
                cells2 = line2.strip().split('\t')
    
                # 按单元格合并数据，这里假设每行的单元格数量相同
                merged_cells = [f"{cell1}:{cell2}" for cell1, cell2 in zip(cells1, cells2)]
    
                # 使用制表符连接单元格，形成新的行
                merged_line = "\t".join(merged_cells)
                merged_lines.append(merged_line)
    
            # 将合并后的数据写入新文件
            with open(output_path, 'w') as outfile:
                for line in merged_lines:
                    outfile.write(line + '\n')
    
        except FileNotFoundError as e:
            print(f"文件未找到: {e}")
        except IOError as e:
            print(f"文件读写错误: {e}")
        except Exception as e:
            print(f"发生错误: {e}")
    
    # 调用函数，传入文件路径和输出文件路径
    merge_files('converted_genotypes-concordance.txt', 'extracted_values-imputed-concordance.txt', 'merged_output-concordance.txt')
    
    def create_matrix_from_line(line):
        # 初始化4x4矩阵
        matrix = [[0.0 for _ in range(4)] for _ in range(4)]
    
        # 移除行尾的换行符并按制表符分割
        entries = line.strip().split('\t')
        
        for entry in entries:
            # 分离索引号和数值字符串
            col_index_str, values_str = entry.split(':')
            col_index = int(col_index_str) 
            values = list(map(float, values_str.split(',')))
    
            # 将数值累加到矩阵中对应的列
            for i, value in enumerate(values):
                if i < 4:  # 确保行索引在4x4矩阵范围内
                    matrix[i][col_index] += value
    
        # 计算总计行和列
        for i in range(3):
            matrix[i][3] = sum(matrix[i][:3])
            matrix[3][i] = sum(matrix[j][i] for j in range(3))
        matrix[3][3] = sum(matrix[i][3] for i in range(3))  # 总计
    
        return matrix
    
    def calculate_Po(matrix):
        # 计算对角线元素的总和
        diagonal_sum = sum(matrix[i][i] for i in range(3))
        # 计算Po
        po = diagonal_sum / matrix[3][3] 
    
        return po  # 确保返回po值
    
    def process_file(file_path, output_file_path):
        with open(file_path, 'r') as file, open(output_file_path, 'w') as output_file:
            lines = file.readlines()
            for line in lines:
                matrix = create_matrix_from_line(line)  # 为每行创建和累加矩阵
                
                # 计算矩阵的po值
                po = calculate_Po(matrix)
                # 打印po值，每个po值一行
                output_file.write(f'{po:.6f}\n')
    
    # 替换为你的文件路径
    input_file_path = 'merged_output-concordance.txt'
    output_file_path = 'concordance_values.txt'  # po值将被保存在这个文件中
    process_file(input_file_path, output_file_path)
    
        
    ##merge chr-pos-R2
    def merge_two_files(file_path1, file_path2, output_file_path):
        with open(file_path1, 'r') as infile1, open(file_path2, 'r') as infile2, open(output_file_path, 'w') as outfile:
            # 读取第一个文件的所有行
            lines1 = infile1.readlines()
            # 读取第二个文件的所有行
            lines2 = infile2.readlines()
            
            # 确保两个文件的行数相同
            if len(lines1) != len(lines2):
                raise ValueError("两个文件的行数不匹配。")
    
            # 合并数据
            for line1, line2 in zip(lines1, lines2):
                columns1 = line1.strip().split('\t')
                columns2 = line2.strip().split('\t')
                
                # 假设每行至少有1列数据
                if len(columns1) >= 1 and len(columns2) >= 1:
                    merged_line = f"{columns1[0]}\t{columns1[1]}\t{columns2[0]}\n"
                    outfile.write(merged_line)
    
    # 使用示例
    file_path1 = 'extracted_chr-pos-concordance.txt'  # 第一个文件的路径
    file_path2 = 'concordance_values.txt'  # 第二个文件的路径
    output_file_path = f'impute5-INFO-{min_threshold:.2f}-{max_threshold:.2f}-chr-pos-concordance-GP-map.txt'  # 输出文件的路径
    merge_two_files(file_path1, file_path2, output_file_path)
