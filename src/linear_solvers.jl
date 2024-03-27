
function GaussElimination(A, b)
    """
        Simple code to perform Gauss elimination consistent with class.

    DOCSTRING

    # Arguments:
    - `A`: DESCRIPTION
    - `b`: DESCRIPTION
    """
    n = length(b)
    Ab = [A b]

    for k = 1:n
        pivot, pivot_row = findmax(abs.(Ab[k:n, k]))
        pivot_row += k - 1
        Ab[[k, pivot_row], :] = Ab[[pivot_row, k], :]

        for i = k+1:n
            factor = Ab[i, k] / Ab[k, k]
            Ab[i, k:n+1] -= factor * Ab[k, k:n+1]
        end
    end

    x = zeros(n)

    for i = n:-1:1
        x[i] = (Ab[i, n+1] - dot(Ab[i, i+1:n], x[i+1:n])) / Ab[i, i]
    end

    return x
end
