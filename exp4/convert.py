def convert_to_markdown_table(input_file_path, output_file_path):
    # 打开输入文件
    with open(input_file_path, "r") as file:
        lines = file.readlines()

    # 提取数据并创建 Markdown 表格
    table_data = []
    for line in lines:
        # 分割行
        parts = line.strip().split()
        if len(parts) == 6:
            method, block_size_x, block_size_y, _, time, _ = parts
            if method == "naive":
                table_data.append([block_size_x, block_size_y, time, ""])
            elif method == "shared_memory":
                # 查找对应的 naive 时间
                for row in table_data:
                    if row[0] == block_size_x and row[1] == block_size_y:
                        row[3] = time
                        break
                else:
                    # 如果没有找到对应的 naive 时间，则添加新行
                    table_data.append([block_size_x, block_size_y, "", time])

    # 创建 Markdown 表格格式
    table = "| block_size_x | block_size_y | naive_time | shared_memory_time |\n"
    table += "| -----------: | -----------: | ---------: | -----------------: |\n"
    for row in table_data:
        table += "| " + " | ".join(row) + " |\n"

    # 将表格写入输出文件
    with open(output_file_path, "w") as file:
        file.write(table)


# 定义文件路径
input_file_path = "./output.txt"
output_file_path = "./table.md"

# 转换并保存 Markdown 表格
convert_to_markdown_table(input_file_path, output_file_path)
