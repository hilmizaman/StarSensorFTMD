def sum_except_index_with_weights(arr, weights, index):
    # Check if the index is valid
    if index < 0 or index >= len(arr):
        return "Invalid index"
    
    # Sum all elements except the one at the given index
    result = sum(val * weight for i, (val, weight) in enumerate(zip(arr, weights)) if i != index)
    
    # Add the missing weight if the index is not at the end of the array
    if index < len(arr) - 1:
        result += arr[index] * weights[index]

    return result

# Example usage
arr = [3, 2, 5, 1, 5, 9]
weights = [1, 1, 1, 1, 1]  # Weights corresponding to each element in arr except the last one
index = 0  # Index of the element to exclude
result = sum_except_index_with_weights(arr, weights, index)
print("Result:", result)
