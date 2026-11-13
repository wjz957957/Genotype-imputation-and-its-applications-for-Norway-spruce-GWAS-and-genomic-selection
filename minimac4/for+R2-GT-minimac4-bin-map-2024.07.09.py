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
        f"INFO/R2>={min_threshold}&INFO/R2<{max_threshold}",
        'minimac4_output_map.vcf.gz'
    ]

    
    # 定义bgzip命令
    bgzip_command = ['bgzip', '-c']
    
    # 将bcftools的输出重定向到bgzip
    # 使用subprocess.Popen()创建bcftools进程，并将输出设置为管道
    bcftools_process = subprocess.Popen(bcftools_view_command, stdout=subprocess.PIPE, text=True)
    
    # 将bcftools的输出作为bgzip的输入
    with open('minimac4_output_R2-0.1-R2.vcf.gz', 'w') as output_file:
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
    input_vcf_gz_path = 'minimac4_output_R2-0.1-R2.vcf.gz'  # 你的.vcf.gz文件路径
    output_file_path = 'extracted_chr-pos-R2.txt'    # 你希望输出的文件路径
    extract_columns_vcf(input_vcf_gz_path, output_file_path)
    
    ##根据提出的染色体位置信息，将重测序里相同的位点提取出来
    
    # 定义bcftools view命令和参数
    bcftools_command = [
        'bcftools', 'view', '--regions-file', 'extracted_chr-pos-R2.txt',
        'Chr11_map_phased_50_rename.vcf.gz'
    ]
    
    # 定义bgzip命令
    bgzip_command = ['bgzip', '-c']
    
    # 创建bcftools进程，并将输出重定向到bgzip的输入
    with open('Chr11_map_phased_50_rename_R2-0.1-R2.vcf.gz', 'wb') as output_file:
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



    ##GTDS-to-DS
    # 指定VCF文件路径
    vcf_filename = 'minimac4_output_R2-0.1-R2.vcf.gz'  # 替换为您的文件名
    
    # 用于存储提取的数据
    extracted_data = []
    
    # 使用 gzip 打开文件并读取内容
    with gzip.open(vcf_filename, 'rt') as vcf_file:
        for line in vcf_file:
            # 跳过注释行
            if line.startswith('##'):
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
    output_filename = 'extracted_columns-imputed-R2.txt'  # 您可以自定义输出文件名
    with open(output_filename, 'w') as output_file:
        for data in extracted_data:
            output_file.write(data + '\n')
    
    # 打印完成消息
    print("Extraction complete. Results are saved in", output_filename)
    
    # 指定输入和输出文件路径
    input_filename = 'extracted_columns-imputed-R2.txt'
    output_filename = 'extracted_values-imputed-R2.txt'
    
    # 读取输入文件并写入输出文件
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        # 逐行读取文件
        for line_number, line in enumerate(infile):
            # 保留第一行
            if line_number == 0:
                outfile.write(line)  # 写入第一行
                continue
    
            # 处理第二行及以后的行
            fields = line.strip().split('\t')  # 根据制表符分割行
            extracted_values = [field.split(':')[0] for field in fields if ':' in field]  # 提取冒号后的值
            extracted_line = '\t'.join(extracted_values)  # 将提取的值重新连接成一行
            outfile.write(extracted_line + '\n')  # 写入结果到输出文件
    
    # 打印完成消息
    print("Extraction complete. Results are saved in", output_filename)
    
    # 指定输入和输出文件路径
    input_filename = 'extracted_values-imputed-R2.txt'  # 原始文件名
    output_filename = 'genotype_012_format-imputed-R2.txt'  # 结果文件名
    
    # 定义一个函数来转换基因型数据
    def convert_genotype(genotype):
        # 将基因型 '0|0' 转换为 '0'，'1|1' 转换为 '2'，'1|0' 和 '0|1' 转换为 '1'
        alleles = genotype.split('|')
        if alleles[0] == '0' and alleles[1] == '0':
            return '0'
        elif alleles[0] == '1' and alleles[1] == '1':
            return '2'
        else:
            return '1'
    
    # 读取输入文件，转换基因型数据，并写入输出文件
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        # 逐行读取输入文件
        for line_number, line in enumerate(infile):
            # 对于第一行标题，直接写入到输出文件
            if line_number == 0:
                outfile.write(line)
            else:
                # 对于数据行，分割每列然后转换基因型
                samples = line.strip().split('\t')
                # 转换每列的基因型数据
                converted_samples = [convert_genotype(sample) for sample in samples]
                # 将转换后的数据写入输出文件，用制表符分隔
                outfile.write('\t'.join(converted_samples) + '\n')
    
    # 打印完成消息
    print("Conversion complete. Results are saved in", output_filename)
    
    ##GT-to-DS
    # 指定VCF文件路径
    vcf_filename = 'Chr11_map_phased_50_rename_R2-0.1-R2.vcf.gz'  # 替换为您的文件名
    
    # 用于存储提取的数据
    extracted_data = []
    
    # 使用 gzip 打开文件并读取内容
    with gzip.open(vcf_filename, 'rt') as vcf_file:
        for line in vcf_file:
            # 跳过注释行
            if line.startswith('##'):
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
    output_filename = 'extracted_columns-R2.txt'  # 您可以自定义输出文件名
    with open(output_filename, 'w') as output_file:
        for data in extracted_data:
            output_file.write(data + '\n')
    
    # 打印完成消息
    print("Extraction complete. Results are saved in", output_filename)
    
    # 指定输入和输出文件路径
    input_filename = 'extracted_columns-R2.txt'  # 原始文件名
    output_filename = 'genotype_012_format-R2.txt'  # 结果文件名
    
    # 定义一个函数来转换基因型数据
    def convert_genotype(genotype):
        # 将基因型 '0|0' 转换为 '0'，'1|1' 转换为 '2'，'1|0' 和 '0|1' 转换为 '1'
        alleles = genotype.split('|')
        if alleles[0] == '0' and alleles[1] == '0':
            return '0'
        elif alleles[0] == '1' and alleles[1] == '1':
            return '2'
        else:
            return '1'
    
    # 读取输入文件，转换基因型数据，并写入输出文件
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        # 逐行读取输入文件
        for line_number, line in enumerate(infile):
            # 对于第一行标题，直接写入到输出文件
            if line_number == 0:
                outfile.write(line)
            else:
                # 对于数据行，分割每列然后转换基因型
                samples = line.strip().split('\t')
                # 转换每列的基因型数据
                converted_samples = [convert_genotype(sample) for sample in samples]
                # 将转换后的数据写入输出文件，用制表符分隔
                outfile.write('\t'.join(converted_samples) + '\n')
    
    # 打印完成消息
    print("Conversion complete. Results are saved in", output_filename)
    
    ##calculate-R2
    def read_data(file_path):
        """读取数据文件，跳过第一行样本名，返回数据列表"""
        data = []
        with open(file_path, 'r') as file:
            next(file, 1)  # 跳过第一行
            for line in file:
                values = line.strip().split('\t')
                data.append([float(value) for value in values[1:]])  # 假设第一列是样本名
        return data
    
    def calculate_pearson_r2(x, y):
            """计算皮尔逊相关系数的平方"""
            if len(x) == 0:
                return np.nan  # 如果没有数据点，返回NaN
            # 计算平均值
            x_bar = np.mean(x)
            y_bar = np.mean(y)
            # 计算总平方和 SSyy
            SSyy = np.sum((y - y_bar)**2)
            
            # 计算协方差 SSxy
            SSxy = np.sum((x - x_bar) * (y - y_bar))
            
            # 计算自变量平方和 SSxx
            SSxx = np.sum((x - x_bar)**2)
            # 确保SSxx不为零，避免除以零
            eps = 1e-15
            SSxx = SSxx + eps
            # 计算回归系数 b_nn
            b_nn = SSxy / SSxx
            
            # 计算残差平方和 SSE
            SSE = SSyy - b_nn * SSxy
            
            # 计算平方相关系数 R_squared
            return  1 - SSE / (SSyy + eps)
    
    def main(file1_path, file2_path, output_path):
        # 读取数据
        data1 = read_data(file1_path)
        data2 = read_data(file2_path)
        
        # 确保两个文件的数据行数相同
        if len(data1) != len(data2):
            raise ValueError("两个文件的数据行数必须相同")
        
        # 计算每行的相关系数平方
        r2_values = []
        for x, y in zip(data1, data2):
            r2 = calculate_pearson_r2(x, y)
            r2_values.append(r2)
        
        # 将结果写入文件
        with open(output_path, 'w') as output_file:
            for r2 in r2_values:
                output_file.write(f"{r2:.4f}\n")  # 格式化输出，保留四位小数
    
    # 使用示例
    file1_path = 'genotype_012_format-imputed-R2.txt'
    file2_path = 'genotype_012_format-R2.txt'
    output_path = 'R2-value.txt'
    
    main(file1_path, file2_path, output_path)
    
    
        
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
    file_path1 = 'extracted_chr-pos-R2.txt'  # 第一个文件的路径
    file_path2 = 'R2-value.txt'  # 第二个文件的路径
    output_file_path = f'minimac4-R2-{min_threshold:.2f}-{max_threshold:.2f}-chr-pos-R2-GT-map.txt'  # 输出文件的路径
    merge_two_files(file_path1, file_path2, output_file_path)