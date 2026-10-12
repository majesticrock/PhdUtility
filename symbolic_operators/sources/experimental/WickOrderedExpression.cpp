#include <mrock/symbolic_operators/experimental/WickOrderedExpression.hpp>
#include <ostream>
#include <algorithm>

namespace mrock::symbolic_operators::experimental {

WickOrderedExpression& WickOrderedExpression::hermitian_conjugate_inplace() noexcept
{
    std::reverse(operators.begin(), operators.end());
    for (auto& op : operators) {
        op.is_daggered = !op.is_daggered;
    }
    return *this;
}

WickOrderedExpression WickOrderedExpression::hermitian_conjugate() const noexcept
{
    WickOrderedExpression copy(*this);
    copy.hermitian_conjugate_inplace();
    return copy;
}


std::ostream& operator<<(std::ostream& os, const WickOrderedExpression& expr) {
    if (expr.empty()) {
        os << ": \\hat{1} :";
    }
    else {
        os << ":" << expr.operators << ":";
    }
    return os;
}

bool operator==(const WickOrderedExpression& lhs, const WickOrderedExpression& rhs) {
    return lhs.operators == rhs.operators;
}

bool operator!=(const WickOrderedExpression& lhs, const WickOrderedExpression& rhs) {
    return !(lhs == rhs);
}
}