while IFS= read -r file_name; do
  # 检查文件名非空
  [ -z "$file_name" ] && continue
  # 执行POA命令并重定向输出
  base_name="${file_name%.*}"  
  /usr/bin/time -v -o "nresult/${base_name}.log" ../POA/bin/POA -i "../POA/data/mt/$file_name" -S > "nresult/${base_name}.msa"
  echo "已处理: $file_name"
done < ids.txt