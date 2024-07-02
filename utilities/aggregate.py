def aggregate(df, group_columns, aggregation_functions): 
  groups = df.groupby(group_columns)
  aggregated = groups.agg(aggregation_functions)  
  df = aggregated.reset_index()
  df.columns = [' '.join(col[::-1]).strip() for col in df.columns.values]
  return df

def aggregate_polars(df, group_columns, aggregation_functions): 
    df_grouped = ( 
        df
        .group_by(group_columns, maintain_order=True)
        .agg(*aggregation_functions)
    )
    return df_grouped
  
