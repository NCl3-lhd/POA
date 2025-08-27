while IFS= read -r file_name; do
  # 检查文件名非空
  [ -z "$file_name" ] && continue
  # 执行POA命令并重定向输出
  ../POA/bin/POA -i "../POA/data/in/$file_name" -p > "./result/bali2dna_test_p/${file_name}.msa"
  echo "已处理: $file_name"
done < ids.txt